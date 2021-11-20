import json
import urllib.request as request
import warnings
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd

from chemlinks import get_chembl, get_uniprot_acc
from drugable.utils import fetch_from_json
from Source import Source

from tqdm import tqdm

class ttdb(Source):
    def __init__(
        self,
        sdir = Path('ttdb/'),
        association_url = 'http://db.idrblab.net/ttd/sites/default/files/ttd_database/P1-07-Drug-TargetMapping.xlsx',
        alias_url = 'http://db.idrblab.net/ttd/sites/default/files/ttd_database/P1-04-Drug_synonyms.txt',
        crossmatch_url = 'http://db.idrblab.net/ttd/sites/default/files/ttd_database/P1-03-TTD_crossmatching.txt',
        up_crossmatch_url = 'http://db.idrblab.net/ttd/sites/default/files/ttd_database/P1-01-TTD_target_download.txt'
        ) -> None:
        super().__init__(
            sdir = sdir,
            association_url = association_url,
            alias_url = alias_url,
            crossmatch_url = crossmatch_url,
            up_crossmatch_url = up_crossmatch_url
            )
        self.drug_fn = self.sdir / 'TTDB_drugs_unified.json'
        self.target_fn = self.sdir / 'TTDB_targets_unified.json'
        self.association_fn = self.sdir / 'TTDB_associations_unified.json'
        self.aliases_fn_raw = self.sdir / 'TTDB_aliases.txt'
        self.aliases_fn_proc = self.sdir / 'TTDB_aliases_proc.csv'
        self.crossmatch_fn_raw = self.sdir / 'TTDB_crossmatching.txt'
        self.crossmatch_fn_proc = self.sdir / 'TTDB_crossmatch_proc.csv'
        self.up_cross_fn_raw = self.sdir / 'TTDB_up_crossmatching.txt'
        self.association_fn_raw = self.sdir / 'TTDB_drug_target_mapping.xlsx'
        self.association_fn_proc = self.sdir / 'TTDB_drug_target_mapping_proc.csv'

    def check_files(self):
        if self.drug_fn.exists() and self.target_fn.exists() and self.association_fn.exists():
            return 0  # fully processed
        elif self.aliases_fn_raw.exists() and self.crossmatch_fn_raw.exists() and self.association_fn_raw.exists():
            return 1  # raw files exist
        return 2  # no files found

    def download(self, association_url = None, aliases_url = None, drug_crossmatch_url = None, up_cross_url = None):
        ttdb_drug_targ_url = association_url or self.association_url
        ttdb_aliases_url = aliases_url or self.alias_url
        ttdb_drug_crossmatch_url = drug_crossmatch_url or self.crossmatch_url
        up_cross_url = up_cross_url or self.up_crossmatch_url

        request.urlretrieve(ttdb_drug_targ_url, self.association_fn_raw)
        request.urlretrieve(ttdb_aliases_url, self.aliases_fn_raw)
        request.urlretrieve(ttdb_drug_crossmatch_url, self.crossmatch_fn_raw)
        request.urlretrieve(up_cross_url, self.up_cross_fn_raw)

    def parse_drugs(self, aliases_fn = None, crossmatch_fn = None, alias_sep='\t', cross_sep='\t'):
        tqdm.pandas(bar_format='{n_fmt} / {total_fmt}')
        aliases_fn = aliases_fn or self.aliases_fn_raw
        crossmatch_fn = crossmatch_fn or self.crossmatch_fn_raw

        crossmatch = pd.read_csv(
            crossmatch_fn,
            skiprows=28,
            sep=cross_sep,
            header=0,
            names=['ttdb_id','id_type','id_value']
            )
        crossmatch = crossmatch.progress_apply(lambda x: x.str.upper().str.strip() if isinstance(x, object) else x)
        crossmatch.to_csv(self.crossmatch_fn_proc)
        aliases = pd.read_csv(
            aliases_fn,
            skiprows=21,
            sep=alias_sep,
            header=0,
            names=['ttdb_id','alias_type','alias_value']
            )
        aliases['alias_value'] = aliases['alias_value'].astype(str)
        aliases = aliases.progress_apply(lambda x: x.str.upper().str.strip() if isinstance(x, object) else x)
        aliases.to_csv(self.aliases_fn_proc)
        cross_gb = crossmatch.groupby(
            ['ttdb_id','id_type'],
            sort=False,
            as_index=False
            )
        cross_agg = cross_gb['id_value'].agg({'id_value' : '|'.join})
        cross_pivot = cross_agg.pivot(
            index='ttdb_id',
            columns='id_type',
            values='id_value'
            )
        cross_pivot.drop(
            ['D_FOMULA','SUPDRCAS','SUPDRATC','TTDDRUID'],
            axis=1,
            inplace=True
            )
        cross_pivot['PUBCHCID'] = cross_pivot['PUBCHCID'].astype(str).progress_apply(lambda x: x.split('; ')[0])
        cross_pivot['PUBCHCID'] = cross_pivot['PUBCHCID'].progress_apply(lambda x: np.nan if x == 'nan' else x)
        alias_gb = aliases.groupby(
            ['ttdb_id',
            'alias_type'],
            sort=False,
            as_index=False
            )
        alias_agg = alias_gb['alias_value'].agg(lambda x: x.unique().tolist())
        alias_pivot = alias_agg.pivot(
            index='ttdb_id',
            columns='alias_type',
            values='alias_value'
            )
        alias_pivot.drop(['DRUGNAME'], axis=1, inplace=True)
        drug_merge = alias_pivot.merge(
            cross_pivot,
            left_on='ttdb_id',
            right_on='ttdb_id'
            )
        drug_merge.rename(
            {
                'CASNUMBE' : 'cas_number',
                'CHEBI_ID' : 'chebi_id',
                'DRUGNAME' : 'ttdb_name',
                'PUBCHCID' : 'pubchem_cid',
                'PUBCHSID' : 'pubchem_sids',
                'TTDDRUID' : 'ttdb_id',
                'SYNONYMS' : 'aliases'
            }, inplace=True, axis=1
        )
        drug_merge['pubchem_sids'] = drug_merge['pubchem_sids'].str.split('; ')
        chembl_dict = get_chembl(drug_merge['pubchem_cid'].tolist(), id_from='cid')
        drug_merge['chembl_id'] = drug_merge['pubchem_cid'].progress_apply(lambda x: chembl_dict[x])
        drug_merge = drug_merge.where(pd.notnull(drug_merge), None)
        drug_merge['ttdb_id'] = drug_merge['ttdb_id'].progress_apply(lambda x: x[0] if x is not None else x)
        fmt_dict_drug = drug_merge.loc[(~drug_merge['chembl_id'].isnull()) | (~drug_merge['pubchem_cid'].isnull())].to_dict('records')
        drug_merge.loc[(drug_merge['chembl_id'].isnull()) &\
                    (drug_merge['pubchem_cid'].isnull())].\
            to_json((self.sdir / 'TTDB_drugs_ununified.json'), 'records')

        with open(self.drug_fn, 'w') as jf:
            d = json.dumps(fmt_dict_drug)
            jf.writelines(d)

    def parse_targets(self, fn = None):
        tqdm.pandas(bar_format='{n_fmt} / {total_fmt}')
        fn = fn or self.up_cross_fn_raw
        up_cross_proc_path = self.sdir / 'TTDB_crossmatch_up_proc.csv'

        with open (fn, 'r') as tf:
            data = tf.readlines()
        data[39] = '\t'.join(['ttdb_id','id_type','id_value','pmid','max_phase']) + '\n'
        with open(fn, 'w') as tf:
            tf.writelines(data)
        up_cross = pd.read_csv(
            fn,
            skiprows=39,
            sep='\t'
            )
        up_cross = up_cross.progress_apply(lambda x: x.str.upper().str.strip() if isinstance(x, object) else x)
        up_cross.to_csv(up_cross_proc_path)

        up_gb = up_cross.groupby(['ttdb_id', 'id_type'], sort=False, as_index=False)
        up_agg = up_gb['id_value'].aggregate(lambda x: x.unique().tolist())
        up_pivot = pd.pivot(up_agg, values='id_value',columns='id_type', index='ttdb_id')
        up_filt = up_pivot[['SYNONYMS','TARGETID','TARGNAME','TARGTYPE','UNIPROID']]
        up_filt.rename(
            {
                'SYNONYMS' : 'aliases',
                'TARGETID' : 'ttdb_id',
                'TARGNAME' : 'ttdb_name',
                'TARGTYPE' : 'ttdb_type',
                'UNIPROID' : 'uniprot_accession'
            }, inplace=True, axis=1
        )
        up_filt['aliases'] = up_filt['aliases'].progress_apply(lambda x: x[0].split('; ') if isinstance(x, list) else '')
        up_filt['ttdb_id'] = up_filt['ttdb_id'].progress_apply(lambda x: x[0])
        up_filt['ttdb_type'] = up_filt['ttdb_type'].progress_apply(lambda x: x[0])
        up_filt['ttdb_name'] = up_filt['ttdb_name'].progress_apply(lambda x: x[0])
        up_filt['uniprot_accession'] =  [[''] if x is np.NaN else x for x in up_filt['uniprot_accession']]
        # TTDB uses the entry_name (XXXX_human) not the entry_id ([P/O]XXXXXX), so must convert
        uniprot_dict = get_uniprot_acc(up_filt.loc[~up_filt['uniprot_accession'].isnull()]['uniprot_accession'].tolist())
        up_filt['uniprot_accession'] = up_filt['uniprot_accession'].progress_apply(lambda x: uniprot_dict[x[0]] if x[0] in uniprot_dict.keys() else '')
        # uniprot link to chembl
        chembl_dict = get_chembl(up_filt['uniprot_accession'].tolist(), id_from='uniprot')
        up_filt['chembl_id'] = up_filt['uniprot_accession'].progress_apply(lambda x: chembl_dict[x] if x in chembl_dict.keys() else '')
        up_filt = up_filt.where(pd.notnull(up_filt), None)
        fmt_dict_up = up_filt.to_dict('records')

        with open(self.target_fn, 'w') as jf:
            d = json.dumps(fmt_dict_up)
            jf.writelines(d)

    def parse_associations(self, fn = None):
        tqdm.pandas(bar_format='{n_fmt} / {total_fmt}')
        fn = fn or self.association_fn_raw
        dt_map = pd.read_excel(fn)
        dt_map = dt_map.progress_apply(lambda x: x.str.upper().str.strip() if isinstance(x, object) else x)
        dt_map.rename(
            {
                'TargetID' : 'ttdb_target_id',
                'DrugID' : 'ttdb_drug_id',
                'MOA' : 'mechanism',
                'Highest_status' : 'max_phase',
                'Activity' : 'activity',
                'Referecnce' : 'reference', # the typo is intentional
            }, inplace=True
        )
        dt_map.to_csv(self.association_fn_proc)

        jdat = dt_map.to_dict(orient='records')
        with open(self.association_fn, 'w') as jf:
            d = json.dumps(jdat)
            jf.writelines(d)

    def load_drugs(self, path=None):
        path = path or self.drug_fn
        return fetch_from_json(path)

    def load_targets(self, path=None):
        path = path or self.target_fn
        return fetch_from_json(path)

    def load_assocations(self, path=None):
        path = path or self.association_fn
        return fetch_from_json(path)

    def run(self, check_files=True):
        print('Running ttdb processing pipeline...')
        if check_files:
            if self.check_files() == 0:
                print("TTDB files exist, skipping")
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
        print('Parsing associations...')
        self.parse_associations()
        return
