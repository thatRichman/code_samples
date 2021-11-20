import gzip
import os
import json
import urllib.request as request
from ast import literal_eval
from pathlib import Path
import ftplib
import ast

import pandas as pd
from drugable.utils.fetch_from_json import fetch_from_json
from drugable.utils.flatten import flatten
from tqdm import tqdm

from Source import Source


class opentargets(Source):
    def __init__(
        self,
        sdir = Path('opentargets/'),
        tract_url = 'http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/input/annotation-files/tractability_buckets-2021-06-03.tsv',
        target_url = 'ftp.ebi.ac.uk/pub/databases/opentargets/platform/21.06/output/etl/json/targets',
        adverse_url = 'http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/input/annotation-files/known_target_safety-2021-06-18.json'
        ) -> None:
        super().__init__(sdir = sdir, tract_url = tract_url, target_url = target_url, adverse_url = adverse_url)
        self.drug_fn = self.sdir / 'opentargets_drugs_unified.json'
        self.target_fn = self.sdir / 'opentargets_targets_unified.json'
        self.tract_fn_raw = self.sdir / 'tractability_file.tsv'
        self.target_fn_raw = self.sdir / 'target_list.csv'
        self.adverse_fn_raw = self.sdir / 'known_target_safety.json'
        self.target_fn_proc = self.sdir / 'target_list.csv'

    def check_files(self):
        if self.drug_fn.exists() and self.target_fn.exists():
            return 0
        elif self.tract_fn_raw.exists() and self.target_fn_raw.exists() and self.adverse_fn_raw.exists():
            return 1
        return 2

    def download(self, tractability_url = None, target_url = None, adverse_url = None):
        # targets are stored on FTP server in multiple JSON
        target_path = Path(target_url) if target_url else Path(self.target_url)
        FTP_URL = target_path.parts[0]
        FTP_USER="anonymous"
        FTP_PASS=""
        ftp = ftplib.FTP(FTP_URL, FTP_USER, FTP_PASS)
        ftp.encoding="utf-8" 
        ftp.cwd(str(target_path.relative_to(target_path.parts[0])))

        for part in ftp.nlst():
            outfile_name = Path(self.sdir) / part
            if not outfile_name.exists():
                outfile = open(outfile_name ,'wb')
                ftp.retrbinary("RETR " + part, outfile.write)
                outfile.close()
        ftp.quit()

        dfs = []
        for part in self.sdir.glob("*.json"):
            dfs.append(pd.read_json(part, lines=True))
            os.remove(part)
        target_df = pd.concat(dfs)
        target_df.to_csv(self.sdir / "target_list.csv", index=False)

        tractability_url = tractability_url or self.tract_url
        adverse_url = adverse_url or self.adverse_url

        request.urlretrieve(tractability_url, self.tract_fn_raw)
        request.urlretrieve(adverse_url, self.adverse_fn_raw)

        # uncompress the target csv
        with gzip.open(self.target_fn_raw, 'rt', newline='') as f:
            with open(self.target_fn_proc, 'wt', newline='') as of:
                of.writelines(f)

    def parse_drugs(self, fn = None, sep='\t'):
        fn = fn or self.tract_fn_raw

        tract_csv = pd.read_csv(fn, sep=sep)
        tract_csv = tract_csv[['symbol','drug_names_dict_sm','drug_names_dict_ab']]
        tract_csv = tract_csv.apply(lambda x: x.str.upper().str.strip() if isinstance(x, object) else x)
        sm_drug_list = [literal_eval(d) for d in tract_csv["drug_names_dict_sm"].tolist()]
        sm_drug_list = [d for d in sm_drug_list if '' not in d.keys()]

        ab_drug_list = [literal_eval(d) for d in tract_csv['drug_names_dict_ab']]
        ab_drug_list = [d for d in ab_drug_list if '' not in d.keys()]

        format_dict = []
        for d in tqdm(sm_drug_list, bar_format='{n_fmt} / {total_fmt}'):
            for k, v in d.items():
                format_dict.append({'chembl_id' : k, 'ot_name' : v})
        for d in tqdm(ab_drug_list, bar_format='{n_fmt} / {total_fmt}'):
            for k, v in d.items():
                format_dict.append({'chembl_id' : k, 'ot_name': v})

        with open(self.drug_fn, 'w') as jf:
            jf.writelines(json.dumps(format_dict))

    def parse_targets(self, targ_fn = None, tract_fn = None, sep=','):
        tqdm.pandas(bar_format='{n_fmt} / {total_fmt}')
        target_fn = targ_fn or self.target_fn_proc
        tract_fn = tract_fn or self.tract_fn_raw

        targ_csv = pd.read_csv(target_fn, sep=sep)
        targ_csv.rename(
            {
                'approvedSymbol' : 'hgnc_symbol',
                'proteinAnnotations' : 'protein_annotations',
                'id' : 'ensembl_id'
            }, axis=1, inplace=True)

        targ_csv["protein_annotations"] = targ_csv["protein_annotations"].apply(lambda x: None if isinstance(x, float) else x)
        targ_csv["protein_annotations"] = targ_csv["protein_annotations"].apply(lambda x: ast.literal_eval(x) if x else None)
        targ_csv["uniprot_accession"] = targ_csv["protein_annotations"].apply(lambda x: x["id"] if x else None)
        targ_csv["add_uniprot_accs"] = targ_csv["protein_annotations"].apply(lambda x: x["accessions"] if x else None)

        targ_csv = targ_csv[['ensembl_id','hgnc_symbol','add_uniprot_accs',"uniprot_accession"]]
        targ_csv = targ_csv.progress_apply(lambda x: x.str.upper().str.strip() if type(x) == str else x)

        tract_csv = pd.read_csv(tract_fn, sep='\t')
        tract_csv.rename(mapper={'ensembl_gene_id' : 'ensembl_id'}, axis=1, inplace=True)
        tract_csv = tract_csv[['ensembl_id', 'accession','symbol']]
        tract_csv['has_tract_dat'] = True
        targ_csv = targ_csv.merge(tract_csv[['symbol','has_tract_dat']], left_on='hgnc_symbol', right_on='symbol', how='outer')
        targ_csv = targ_csv[['ensembl_id','hgnc_symbol','add_uniprot_accs','uniprot_accession','has_tract_dat']]
        targ_csv.to_json(self.target_fn, orient='records')

    def parse_adverse(self, fn = None):
        #TODO
        pass

    def load_drugs(self, path=None):
        path = path or self.drug_fn
        return fetch_from_json(path)

    def load_targets(self, path=None):
        path = path or self.target_fn
        return fetch_from_json(path)

    def run(self, check_files=True):
        print('Running opentargets processing pipeline...')
        if check_files:
            if self.check_files() == 0:
                print('OpenTargets files exist, skipping')
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
        print('Parsing adverse events [TODO]...')
        self.parse_adverse()
        return
