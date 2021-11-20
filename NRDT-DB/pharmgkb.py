import csv
import json
import urllib.request as request
import zipfile
from ast import literal_eval
from pathlib import Path

import numpy as np
import pandas as pd
from drugable.utils.fetch_from_json import fetch_from_json
from drugable.utils.flatten import flatten
from tqdm import tqdm

from chemlinks import get_chembl, get_phase
from Source import Source


class pharmgkb(Source):

    def __init__(
        self,
        sdir = Path('pharmgkb/'),
        drug_url = 'https://api.pharmgkb.org/v1/download/file/data/drugs.zip',
        target_url = 'https://api.pharmgkb.org/v1/download/file/data/genes.zip') -> None:
        super().__init__(sdir = sdir,  drug_url = drug_url, target_url = target_url)
        self.drug_fn = self.sdir / 'pharmgkb_drugs_unified.json'
        self.target_fn = self.sdir / 'pharmgkb_targets_unified.json'
        self.drug_fn_raw = self.sdir / 'pgkb_drugs.zip'
        self.target_fn_raw = self.sdir / 'pgkb_targets.zip'
        self.drug_fn_proc = self.sdir / 'drugs.tsv'
        self.target_fn_proc = self.sdir / 'genes.tsv'

    def check_files(self):
        if self.drug_fn.exists() and self.target_fn.exists():
            return 0
        elif self.drug_fn_raw.exists() and self.target_fn_raw.exists():
            return 1
        return 2

    def _parse_cross_refs(self, cf):
        t = cf.replace('"','').replace('\"','').split(',')
        v = [e.split(':',maxsplit=1) for e in t]
        cross_dict = {}
        for e in v:
            if e[0] in cross_dict.keys():
                cross_dict[e[0]].append(e[1])
            else:
                cross_dict[e[0]] = [e[1]]
        for k, v in cross_dict.items():
            if isinstance(v, list) and len(v) == 1:
                cross_dict[k] = v[0]
        return cross_dict

    def _add_cd(self, row):
        if row['Cross-references'] is not None:
            cd = self._parse_cross_refs(row['Cross-references'])
            for k, v in cd.items():
                row[k] = v
            return row
        else:
            return row

    def download(self, drug_url = None, target_url = None):
        drug_url = drug_url or self.drug_url
        target_url = target_url or self.target_url

        opener = request.build_opener()
        opener.addheaders = [('User-Agent', 'User-AgentMozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:84.0) Gecko/20100101 Firefox/84.0')]
        request.install_opener(opener)
        request.urlretrieve(drug_url, self.drug_fn_raw)
        request.urlretrieve(target_url, self.target_fn_raw)

        with zipfile.ZipFile(self.drug_fn_raw, 'r') as zf:
            zf.extract('drugs.tsv', self.sdir)
        with zipfile.ZipFile(self.target_fn_raw, 'r') as zf:
            zf.extract('genes.tsv', self.sdir)

    def parse_drugs(self, fn = None):
        tqdm.pandas(bar_format='{n_fmt} / {total_fmt}')
        fn = fn or self.drug_fn_proc
        drug_tsv = pd.read_csv(fn, sep='\t', quoting = csv.QUOTE_ALL)
        drug_tsv = drug_tsv.fillna(np.nan).replace([np.nan], [None])
        drug_tsv['max_phase'] = drug_tsv['Top Clinical Annotation Level'].progress_apply(lambda x: get_phase(x))
        drug_tsv['Generic Names'] = drug_tsv['Generic Names'].progress_apply(lambda x: x.replace('"','').split(',') if x else None)
        drug_tsv['Trade Names'] = drug_tsv['Trade Names'].progress_apply(lambda x: x.replace('"','').split(',') if x else None)
        drug_tsv['ATC Identifiers'] = drug_tsv['ATC Identifiers'].progress_apply(lambda x: x.replace('"','').split(',') if x else None)
        drug_tsv['aliases'] = drug_tsv.progress_apply(lambda x: list(set(flatten([x['Name'], x['Generic Names'], x['Trade Names']]))), axis=1)
        # parse cross references
        drug_tsv = drug_tsv.apply(lambda x: self._add_cd(x), axis=1)
        drug_tsv['uniprot_accession'] = drug_tsv['UniProtKB'].progress_apply(lambda x: x[0] if isinstance(x, list) else x)
        drug_tsv = drug_tsv.fillna(np.nan).replace([np.nan], [None])
        drug_tsv.rename(
            {
                'ChEBI' : 'chebi_id',
                'ChemSpider' : 'chemspider_id',
                'BindingDB' : 'bindingdb_id',
                'Chemical Abstracts Service' : 'cas_id',
                'UniProtKB' : 'add_uniprot_accs',
                'Therapeutic Targets Database' :  'ttdb_id',
                'National Drug Code Directory' : 'ndcd_id',
                'RxNorm Identifiers' : 'add_rxnorm_ids',
                'KEGG Drug' : 'kegg_drug_id',
                'KEGG Compound' : 'kegg_compound_id',
                'HMDB' : 'hmdb_id',
                'GenBank' : 'genbank_id',
                'PubChem Substance' : 'pubchem_sids',
                'PharmGKB Accession Id' : 'pharmgkb_id',
                'Name' : 'pharmgkb_name',
                'Generic Names' : 'generic_names',
                'Trade Names' : 'trade_names',
                'Type' : 'pgkb_type',
                'Cross-references' : 'cross_references',
                'ATC Identifiers' : 'atc_ids',
                'PubChem Compound Identifiers' : 'pubchem_cid'
            }, axis = 1, inplace =  True
        )
        drug_tsv['add_chebi_ids'] = drug_tsv['chebi_id']
        drug_tsv['chebi_id'] = drug_tsv['chebi_id'].progress_apply(lambda x: x[0] if isinstance(x, list) else x)
        drug_tsv['add_ttdb_ids'] = drug_tsv['ttdb_id']
        drug_tsv['ttdb_id'] = drug_tsv['ttdb_id'].progress_apply(lambda x: x[0] if isinstance(x, list) else x)
        drug_tsv = drug_tsv[['pharmgkb_id','aliases','pgkb_type','chebi_id','chemspider_id','bindingdb_id','cas_id','add_uniprot_accs',
                            'ttdb_id','ndcd_id','add_rxnorm_ids','kegg_drug_id','kegg_compound_id','hmdb_id','genbank_id','pubchem_sids',
                            'pubchem_cid','atc_ids','max_phase','uniprot_accession','add_chebi_ids']]

        drug_tsv.to_json(self.drug_fn, orient='records')
        drug_tsv = drug_tsv.applymap(lambda x: x.upper() if type(x) == str else x)
        drug_tsv = drug_tsv.applymap(lambda x: [y.upper() for y in x] if type(x) == list else x)

    def parse_targets(self, fn=None):
        tqdm.pandas(bar_format='{n_fmt} / {total_fmt}')
        fn = fn or self.target_fn_proc
        targ_tsv = pd.read_csv(fn, sep='\t', quoting=csv.QUOTE_ALL)
        targ_tsv = targ_tsv.fillna(np.nan).replace([np.nan], [None])
        targ_tsv['Alternate Names'] = targ_tsv['Alternate Names'].progress_apply(lambda x: x.replace('"','').upper().split(',') if x else None)
        targ_tsv['Alternate Symbols'] = targ_tsv['Alternate Symbols'].progress_apply(lambda x: x.replace('"','').upper().split(',') if x else None)
        targ_tsv['aliases'] = targ_tsv.progress_apply(lambda x: list(set(flatten([x['Name'].upper(), x['Alternate Names'],x['Alternate Symbols']]))), axis=1)
        targ_tsv = targ_tsv.progress_apply(lambda x: self._add_cd(x), axis=1)
        targ_tsv = targ_tsv.fillna(np.nan).replace([np.nan], [None])
        targ_tsv['uniprot_accession'] = targ_tsv['UniProtKB'].progress_apply(lambda x: x[0] if isinstance(x, list) else x)
        targ_tsv.rename(
            {
                'ALFRED' : 'alfred_id',
                'NCBI Gene ID' : 'ncbi_gene_id',
                'Ensembl Id' : 'ensembl_id',
                'PharmGKB Accession Id' : 'pharmgkb_id',
                'Symbol' : 'hgnc_symbol',
                'UniProtKB' : 'add_uniprot_accs',
                'RefSeq Protein' : 'refseq_prot_id',
                'RefSeq RNA' : 'refseq_rna_id',
                'RefSeq DNA' : 'refseq_dna_id',
                'HumanCyc Gene' : 'humancyc_gene_id',
                'Comparative Toxicogenomics Database' : 'comp_tox_id',
                'OMIM' : 'omim_id',
                'GenAtlas' : 'gen_atlas_id',
            }, axis = 1, inplace = True
        )
        targ_tsv = targ_tsv[['alfred_id','ncbi_gene_id','ensembl_id','aliases','pharmgkb_id','hgnc_symbol','add_uniprot_accs',
                            'refseq_prot_id','refseq_rna_id','refseq_dna_id','humancyc_gene_id','comp_tox_id','omim_id',
                            'gen_atlas_id','uniprot_accession']]
        targ_tsv = targ_tsv.fillna(np.nan).replace([np.nan], ['NONE'])
        targ_tsv = targ_tsv.progress_apply(lambda x: x.str.upper() if isinstance(x, object) else x)
        targ_tsv = targ_tsv.fillna(np.nan).replace(['NONE'], [None])
        chembl_dict = get_chembl(targ_tsv['uniprot_accession'].tolist(), id_from='uniprot')

        targ_tsv['chembl_id'] = targ_tsv['uniprot_accession'].progress_apply(lambda x: chembl_dict[x] if x in chembl_dict.keys() else None)
        targ_tsv.to_json(self.target_fn, orient='records')

    def load_drugs(self, path = None) -> list:
        if not path:
            path = self.drug_fn
        return fetch_from_json(path)

    def load_targets(self, path = None) -> list:
        if not path:
            path = self.target_fn
        return fetch_from_json(path)

    def run(self, check_files=True):
        print('Running pharmgkb processing pipeline...')
        if check_files:
            if self.check_files() == 0:
                print("PharmGKB files exist, skipping")
                return
            elif self.check_files() == 2:
                print('Downloading data...')
                self.download()
        else:
            print('Downloading data...')
            self.download()
        print('Parsing drugs...')
        self.parse_drugs()
        print('Parsing targets...')
        self.parse_targets()
        return
