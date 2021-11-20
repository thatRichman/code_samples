import datetime
import sqlite3
from pathlib import Path
import pandas as pd

class DrugDB:
    def __init__(self, path=None):
        self.ts = datetime.datetime.today().strftime('%Y%m%d_%H%M%S')
        if not path:
            path = Path(self.ts + '_drug_database.db')
        else:
            path = Path(path)
        self.conn = sqlite3.connect(path)
        self.cursor = self.conn.cursor()
        self.cursor.executescript(
        """
CREATE TABLE IF NOT EXISTS drug (
    uuid BLOB NOT NULL PRIMARY KEY,
    chembl_id TEXT,
    ps_id TEXT,
    pharmgkb_id TEXT,
    ttdb_id TEXT,
    comp_tox_id TEXT,
    pubchem_cid TEXT,
    pubchem_sids TEXT,
    cas_number TEXT,
    chebi_id TEXT,
    aliases TEXT,
    max_phase TEXT,
    chembl_type TEXT,
    pharmgkb_type TEXT,
    ttdb_type TEXT,
    has_tract_dat INT,
    source TEXT
); Jeffrey Lengyel
    chembl_id TEXT,
    ensembl_id TEXT,
    ps_id TEXT,
    ttdb_id TEXT,
    gen_atlas_id TEXT,
    humancyc_gene_id TEXT,
    ncbi_gene_id TEXT,
    omim_id TEXT,
    pharmgkb_id TEXT,
    refseq_dna_id TEXT,
    refseq_prot_id TEXT,
    refseq_rna_id TEXT,
    alfred_id TEXT,
    uniprot_accession TEXT,
    add_uniprot_accs TEXT,
    hgnc_symbol TEXT,
    other_symbols TEXT,
    aliases TEXT,
    organism TEXT,
    chembl_type TEXT,
    ttdb_type TEXT,
    pharmgkb_type TEXT,
    source TEXT
);
CREATE TABLE IF NOT EXISTS association (
    uuid BLOB NOT NULL PRIMARY KEY,
    drug_id BLOB NOT NULL,
    target_id BLOB NOT NULL,
    mechanism INTEGER NOT NULL,
    mechanism_text TEXT,
    mechanism_original,
    direct_interaction INTEGER,
    conflict_resolved INTEGER,
    source TEXT NOT NULL,
    CONSTRAINT drug_id_fk FOREIGN KEY (drug_id) REFERENCES drug(uuid),
    CONSTRAINT target_id_fk FOREIGN KEY (target_id) REFERENCES target(uuid)
);
CREATE TABLE IF NOT EXISTS adverse (
    uuid BLOB NOT NULL PRIMARY KEY,
    target_id TEXT,
    adverse_events JSON,
    CONSTRAINT target_id_fk FOREIGN KEY (target_id) REFERENCES target(rowid)
);
        """
        )
        self.conn.commit()

    def __del__(self):
        self.conn.commit()
        self.conn.close()

    def build(self, drugs, targets, associations):
        self.cursor.executemany(
        """
INSERT INTO drug (
    uuid,
    chembl_id,
    ps_id,
    pharmgkb_id,
    ttdb_id,
    comp_tox_id,
    pubchem_cid,
    pubchem_sids,
    cas_number,
    chebi_id,
    aliases,
    max_phase,
    chembl_type,
    pharmgkb_type,
    ttdb_type,
    has_tract_dat,
    source
) VALUES (
    :uuid,
    :chembl_id,
    :ps_id,
    :pharmgkb_id,
    :ttdb_id,
    :comp_tox_id,
    :pubchem_cid,
    :pubchem_sids,
    :cas_number,
    :chebi_id,
    :aliases,
    :max_phase,
    :chembl_type,
    :pharmgkb_type,
    :ttdb_type,
    :has_tract_dat,
    :src
    )
        """, drugs
        )
        self.cursor.executemany(
        """
INSERT INTO target (
    uuid,
    chembl_id,
    ensembl_id,
    ps_id,
    ttdb_id,
    gen_atlas_id,
    humancyc_gene_id,
    ncbi_gene_id,
    omim_id,
    pharmgkb_id,
    refseq_dna_id,
    refseq_prot_id,
    refseq_rna_id,
    alfred_id,
    uniprot_accession,
    add_uniprot_accs,
    hgnc_symbol,
    other_symbols,
    aliases,
    organism,
    chembl_type,
    pharmgkb_type,
    ttdb_type,
    source
) VALUES (
    :uuid,
    :chembl_id,
    :ensembl_id,
    :ps_id,
    :ttdb_id,
    :gen_atlas_id,
    :humancyc_gene_id,
    :ncbi_gene_id,
    :omim_id,
    :pharmgkb_id,
    :refseq_dna_id,
    :refseq_prot_id,
    :refseq_rna_id,
    :alfred_id,
    :uniprot_accession,
    :add_uniprot_accs,
    :hgnc_symbol,
    :other_symbols,
    :aliases,
    :organism,
    :chembl_type,
    :pharmgkb_type,
    :ttdb_type,
    :src
)
        """, targets
        )
        self.cursor.executemany(
        """
INSERT INTO association (
    uuid,
    drug_id,
    target_id,
    mechanism,
    mechanism_text,
    mechanism_original,
    direct_interaction,
    conflict_resolved,
    source
) VALUES (
    :uuid,
    :drug_uuid,
    :target_uuid,
    :mechanism,
    :mechanism_text,
    :mechanism_original,
    :direct_interaction,
    :conflict_resolved,
    :src
)
        """, associations
        )

    def _construct_filter(k, v):
        if isinstance(v, (list, tuple)):
            return f"AND {k} IN ({','.join(v)})"
        elif isinstance(v, bool):
            if v:
                return f'AND {k} == 1'
            else:
                return f'AND {k} == 0'
        else:
            return f'AND {k} == {v}'

    def select_drugs_by(self, type, values, db_from=None, assoc_filters = None, **kwargs):
        """Build and execute query to return drug-target relations from drug database

        Args:
            type (str): one of 'target', 'class', 'alias', 'id'
            values (iterable): values corresponding to selected type
            db_from (str, optional): database that values are derived from. Defaults to None.
            assoc_filters (dict, optional): optional dictionary of association filters. Defaults to None.

            db_from is required when type == 'id' and recommended when type == 'target'
            If db_from is not specified when type == 'target', will search for target aliases

            type == 'alias' matches by drug aliases, not target aliases

            to filter the inner query (i.e., to filter on columns of either drug or target tables), pass
            filters as keyword arguments e.g. 'max_phase' = 4, or organism = 'Homo sapiens'
        """
        # turn kwarg filters in query chunks
        filters = []
        for k, v in kwargs.items():
            filters.append(self._construct_filter(k,v))
        filt_str = " ".join(filters)

        # optionally filter returned associations
        a_filters = []
        if assoc_filters:
            for k, v in assoc_filters.items():
                a_filters.append(self._construct_filter(k,v))
        a_filt_str = " ".join(a_filters)

        if type == 'target':
            # return all drugs that act on each target (alias) in values
            if not db_from:
                db_from = 'aliases'
            query_base = f"""
            SELECT *
            FROM association
            WHERE target_uuid
            IN (
                SELECT uuid
                FROM target
                WHERE {db_from} in :values
                """
        if type == 'class':
            # return all drugs in each class in values
            query_base = """
            SELECT *
            FROM ASSOCIATION
            WHERE target_uuid
            IN (
                SELECT uuid
                FROM target
                WHERE (chembl_type IN :values
                OR ttdb_type IN :values
                OR pharmgkb_type IN :values)
            """
        if type == 'alias':
            # select drugs with alias in values
            query_base = f"""
            SELECT *
            FROM association
            WHERE drug_uuid
            IN (
                SELECT uuid
                FROM drug
                WHERE aliases in :values
            """
        if type == 'id':
            if not db_from:
                raise ValueError('db_from is required when selecting drugs by id')
            query_base = f"""
            SELECT *
            FROM association
            WHERE drug_uuid
            IN (
                SELECT uuid
                FROM drug
                WHERE {db_from} in :values
                """
        res = self.cursor.execute(query_base + filt_str + ')' + a_filt_str)
        cols = [column[0] for column in res.description]
        return cols, res.fetchall()

    def background_data(self, entities, db_from=None):
        """convienent wrapper around select_drugs_by for returning the total set of background drugs
        for a given set of entities.

        Background data should be as complete as possible, so excessive filtering is discouraged.
        If you must filter the background data, call select_drugs_by directly

        Args:
            targets ([type]): [description]
            db_from ([type], optional): [description]. Defaults to None.
        """
        return self.select_drugs_by(type='target', values = entities, db_from = db_from)