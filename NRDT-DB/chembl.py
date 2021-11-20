from pathlib import Path
import json

from chembl_webresource_client.new_client import new_client

from drugable.utils.fetch_from_json import fetch_from_json
from chemlinks import CHEMBL_DIR, get_phase, get_cid, get_mechanism
from Source import Source
from tqdm import tqdm

class chembl(Source):
    def __init__(
        self,
        sdir = Path('chembl/'),
        ) -> None:
        super().__init__(sdir = sdir)
        self.drug_fn = self.sdir / 'chembl_drugs_unified.json'
        self.target_fn = self.sdir / 'chembl_targets_unified.json'
        self.association_fn = self.sdir / 'chembl_associations_unified.json'

    def check_files(self):
        if self.drug_fn.exists() and self.target_fn.exists() and self.association_fn.exists():
            return 0
        return 2

    def parse_drugs(self):
        """ fetch chembl data"""
        drug = new_client.drug
        drug.set_format('json')
        res_drug = drug.all().only('molecule_chembl_id')
        chembl_ids = [r['molecule_chembl_id'] for r in res_drug]
        molec = new_client.molecule
        res_molec = molec.all().filter(molecule_chembl_id__in=chembl_ids).only(
            'pref_name',
            'molecule_synonyms',
            'molecule_chembl_id',
            'molecule_type',
            'cross_references',
            'max_phase'
            )
        sids = []
        for i in tqdm(res_molec, bar_format='{n_fmt} / {total_fmt}'):
            sids.extend([j['xref_id'] for j in i['cross_references'] if j['xref_src'] == 'PubChem']) # bottleneck, unavoidable
        cid_ref_dict = get_cid(list(set(sids)), id_from='sid')

        format_dict = []
        for i in tqdm(res_molec, bar_format='{n_fmt} / {total_fmt}'):
            pubchem_sids = [j['xref_id'] for j in i['cross_references'] if j['xref_src'] == 'PubChem']
            try:
                pubchem_cid = cid_ref_dict[pubchem_sids[0]]
            except IndexError as ie: # no available SIDs
                pubchem_cid = ''
            except KeyError as ke: # didn't find matching CID (shouldn't really happen)
                pubchem_cid = ''
            format_dict.append({
                'name' : i['pref_name'].upper(),
                'aliases' : list(map(str.upper, set([j['molecule_synonym'] for j in i['molecule_synonyms']]))),
                'chembl_id' : i['molecule_chembl_id'].upper(),
                'chembl_type' : i['molecule_type'].upper(),
                'pubchem_sids' : list(map(str.upper, set(pubchem_sids))),
                'pubchem_cid' : str(pubchem_cid).upper(),
                'max_phase' : get_phase(i['max_phase'])
                })

        with open(self.drug_fn, 'w') as jf:
            j = json.dumps(format_dict)
            jf.writelines(j)

    def parse_targets(self):
        target = new_client.target
        res_targ = target.all().only(
            'pref_name',
            'target_type',
            'organism',
            'target_components',
            'target_chembl_id')
        format_dict = []
        for i in tqdm(res_targ, bar_format='{n_fmt} / {total_fmt}'):
            try:
                aliases = '|'.join([s['component_synonym'] for s in i['target_components'][0]['target_component_synonyms']])
            except (IndexError, KeyError) as ie:
                aliases = ''
            try:
                symbol = [s['component_synonym'] for s in i['target_components'][0]['target_component_synonyms'] if s['syn_type'] == 'GENE_SYMBOL'][0]
            except (IndexError, KeyError) as ie:
                symbol = ''
            try:
                ensembl_ids = [j['xref_id'] for j in i['target_components'][0]['target_component_xrefs'] if j['xref_src_db'] == 'EnsemblGene']
                ensembl_id = ensembl_ids[0]
            except (IndexError, KeyError) as ie:
                ensembl_ids = ['']
                ensembl_id = ''
            try:
                uniprot_accessions = [j['xref_id'] for j in i['target_components'][0]['target_component_xrefs'] if j['xref_src_db'] == 'UniProt']
                uniprot_acc = uniprot_accessions[0]
            except (IndexError, KeyError) as ie:
                uniprot_accessions = ['']
                uniprot_acc = ''

            if str(i['target_type']).upper() not in ['ORGANISM','TISSUE'] and str(i['organism']).upper() == 'HOMO SAPIENS':
                format_dict.append({
                    'aliases' : aliases.upper(),
                    'hgnc_symbol' : symbol.upper(),
                    'chembl_name' : i['pref_name'].upper(),
                    'chembl_type' : i['target_type'].upper(),
                    'organism' : str(i['organism']).upper(),
                    'chembl_id' : i['target_chembl_id'].upper(),
                    'ensembl_id' : ensembl_id.upper(),
                    'add_ensembl_ids' : list(map(str.upper, set(ensembl_ids))),
                    'uniprot_accession' : uniprot_acc.upper(),
                    'add_uniprot_accs' : list(map(str.upper, set(uniprot_accessions)))
                    })

        with open(self.target_fn, 'w') as jf:
            j = json.dumps(format_dict)
            jf.writelines(j)

    def parse_associations(self):
        # parse chembl drug-target interactions
        moa = new_client.mechanism
        moa.set_format('json')
        all_mechs = moa.all()
        format_dict = []
        for x in tqdm(all_mechs, bar_format='{n_fmt} / {total_fmt}'):
            tmp = {}
            tmp['mechanism'], tmp['mechanism_original'] = get_mechanism(x['action_type'])
            tmp['mechanism_text'] = x['mechanism_of_action'].upper() if x['mechanism_of_action'] else x['mechanism_of_action']
            tmp['drug_chembl_id'] = x['molecule_chembl_id'].upper() if x['molecule_chembl_id'] else x['molecule_chembl_id']
            tmp['target_chembl_id'] = x['target_chembl_id'].upper() if x['target_chembl_id'] else x['target_chembl_id']
            tmp['direct_interaction'] = x['direct_interaction']
            format_dict.append(tmp)

        with open(self.association_fn, 'w') as jf:
                j = json.dumps(format_dict)
                jf.writelines(j)

    def load_drugs(self, path = None):
        path = path or self.drug_fn
        return fetch_from_json(path)

    def load_targets(self, path = None):
        path = path or self.target_fn
        return fetch_from_json(path)

    def load_associations(self, path = None):
        path = path or self.association_fn
        return fetch_from_json(path)

    def run(self, check_files=True):
        print('Running chembl processing pipeline...')
        if check_files:
            if self.check_files() == 0:
                print("chembl files exist, skipping.")
                return
        # unprocessed chembl data is not stored locally
        print('Parsing drugs...')
        self.parse_drugs()
        print('Parsing targets...')
        self.parse_targets()
        print('Parsing associations...')
        self.parse_associations()
        return
