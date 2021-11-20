import gzip
import json
import urllib.request as request
import warnings
from io import StringIO
from pathlib import Path
from time import sleep
from urllib import parse

import pandas as pd
from chembl_webresource_client.new_client import new_client
from drugable.utils.flatten import flatten
from drugable.utils.list_chunk import list_chunk

from tqdm import tqdm

from math import ceil

CHEMBL_DIR = Path('chembl/')

PUG_REST_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{}/cids/JSON?cids_type=standardized'
UNICHEM_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz'
UNIPROT_URL = 'https://www.uniprot.org/uploadlists/'

INHIBITION_KEYS = ('INHIBITION','inhibition','INHIBITOR', 'inhibitor', 'ANTAGONIST','antagonist','DOWN REGULATE', 'down regulate', 'DOWN REGULATION', 'down regulation',
                    'PARTIAL ANTAGONIST','partial antagonist', 'negative allosteric modulator', 'NEGATIVE ALLOSTERIC MODULATOR', 'ALLOSTERIC ANTAGONIST',
                    'ANTISENSE INHIBITOR', 'INVERSE AGONIST', 'BLOCKER','ANTISENSE INHIBITOR','RNAI INHIBITOR','NEGATIVE MODULATOR') # map to -1
ACTIVATION_KEYS = ('AGONIST', 'agonist', 'UPREGULATOR','upregulator', 'UP REGULATION', 'up regulation',
                    'UP REGULATE','up regulate','PROMOTER','promoter','ACTIVATOR','activator',
                    'PARTIAL AGONIST','partial agonist', 'POSITIVE ALLOSTERIC MODULATOR',
                    'positive allosteric modulator', 'RNAI ACTIVATOR','INVERSE ANTAGONIST','ANTISENSE ACTIVATOR',
                    'POSITIVE MODULATOR') # map to 1
AMBIGUOUS_KEYS = ('MODULATOR','modulator', 'BINDING AGENT', 'binding agent') # map to 0

PHASE_DICT = {
    'PHASE1' : 'PHASE 1', 'phase1' : 'PHASE 1', 'PHASE 1' : 'PHASE 1', '1' : 'PHASE 1', 1 : 'PHASE 1', 'one' : 'PHASE 1', 'phase 1' : 'PHASE 1', 'Phase 1' : 'PHASE 1',
    'PHASE2' : 'PHASE 2', 'phase2' : 'PHASE 2', 'PHASE 2' : 'PHASE 2', '2' : 'PHASE 2', 2 : 'PHASE 2', 'two' : 'PHASE 2', 'phase 2' : 'PHASE 2', 'Phase 2' : 'PHASE 2',
    'PHASE3' : 'PHASE 3', 'phase3' : 'PHASE 3', 'PHASE 3' : 'PHASE 3', '3' : 'PHASE 3', 3 : 'PHASE 3', 'three' : 'PHASE 3', 'phase 3' : 'PHASE 3', 'Phase 3' : 'PHASE 3',
    'PHASE4' : 'PHASE 4', 'phase4' : 'PHASE 4', 'PHASE 4' : 'PHASE 4', '4' : 'PHASE 4', 4 : 'PHASE 4', 'four' : 'PHASE 4', 'phase 4' : 'PHASE 4', 'Phase 4' :'PHASE 4',
    'investigational' : 'INVESTIGATIONAL', 'Investigational' : 'INVESTIGATIONAL', 'INVESTIGATIONAL' : 'INVESTIGATIONAL', 'Investigative' : 'INVESTIGATIONAL', 'investigative' : 'INVESTIGATIONAL',
    'approved' : 'APPROVED', 'APPROVED' : 'APPROVED', 'Approved' : 'APPROVED', 'Marketed' : 'APPROVED', 'MARKETED' : 'APPROVED', 'indicated' : 'APPROVED', 'INDICATED' : 'APPROVED',
    'Patented' : 'APPROVED', 'Discontinued' : 'DISCONTINUED', 'Discontinued in Phase 1' : 'DISCONTINUED', 'Discontinued in Phase 2' : 'DISCONTINUED', 'Discontinued in Phase 3' : 'DISCONTINUED',
    'terminated' : 'TERMINATED', 'TERMINATED' : 'TERMINATED', 'Terminated' : 'TERMINATED','WITHDRAWN' : 'TERMINATED', 'withdrawn' : 'TERMINATED', 'Withdrawn' : 'TERMINATED',
    'Withdrawn from market' : 'WITHDRAWN FROM MARKET', 'preclinical' : 'PRECLINICAL', 'PRECLINICAL' : 'PRECLINICAL', 'Preclinical' : 'PRECLINICAL',
    '0' : '0', '1A' : 'PHASE 1A', '1B' : 'PHASE 1B', '2A' : 'PHASE 2A','2B' : 'PHASE 2B', None : None, 'None' : None, 'NONE': None
}


def get_mechanism(x):
    """Attempt to map mechanism string onto correct mechanism value (-1, 1, or 0)

    If x is not in keys, returns x with a warning

    Args:
        x (str): mechanism string

    Returns:
        string: mapped mechanism
    """
    if x in INHIBITION_KEYS:
        return -1, x
    elif x in ACTIVATION_KEYS:
        return 1, x
    elif x in AMBIGUOUS_KEYS:
        warnings.warn(f'Ambiguous mechanism of action {x}, assigning 0')
        return 0, x
    else:
        warnings.warn(f'Unrecognized mechanism of action {x}, assigning 0')
        return 0, x


def get_phase(x):
    """Attempt to map phase string onto correct phase representation

    If x cannot be mapped, returns x silently

    Args:
        x (str): phase string to be mapped

    Returns:
        str: mapped phase
    """
    return PHASE_DICT.get(x, x)


def get_chembl_pc_linker(path):
    """Downloads and unzips Unichem's pubchem-->chembl ID linker database

    Args:
        path (str or Path): directory path to place the file
    """
    # retrieve the file
    unichem_path = (Path(path) / 'unichem_linker.json') or (CHEMBL_DIR / 'unichem_linker.json')
    unichem_fn = (Path(path) / 'unichem_linker.txt.gz') or (CHEMBL_DIR / 'unichem_linker.txt.gz')
    request.urlretrieve(UNICHEM_URL, unichem_fn)
    # unzip the file
    with gzip.open(unichem_fn, 'rt', newline='') as f:
        ugz_fn = Path(str(unichem_fn).replace(".gz", ""))
        with open(ugz_fn, 'wt', newline='') as of:
            of.writelines(f)
    # convert two-column tsv to json object
    ucdf = pd.read_csv(ugz_fn, sep='\t', header=0, names=['chembl_id','pubchem_id'])
    ucdf['pubchem_id'] = ucdf['pubchem_id'].astype(str)
    ucdf.to_json(unichem_path, orient='records')
    unichem_fn.unlink() # delete gz version
    ugz_fn.unlink() # delete txt version


def load_chembl_pc_linker(linker_fn):
    """Load the pubchem-->chembl linker into memory

    Args:
        linker_fn (str or Path): path tolinker file

    Returns:
        dict: pubchem_id : chembl_id dict
    """
    if not Path(linker_fn).exists():
        get_chembl_pc_linker(Path(linker_fn).parents[0])
    with open(linker_fn, 'r') as jf:
        j = json.load(jf)
    j_dict = {e['pubchem_id'] : e['chembl_id'] for e in j}
    return j_dict


def get_chembl(x, id_from='cid', linker_fn = None):
    """Convert various IDs to chembl id

    Currently supports pubchem CIDs ('cid') and uniprot IDs ('uniprot')

    Args:
        x (str or list): IDs to convert
        id_from (str, optional): database to convert from. Defaults to 'cid'.
        linker_fn (str or Path, optional): path to pubchem-->chembl linker. Defaults to None. If None, searches in default location

    Raises:
        ValueError: Raised when unspported id_from is provided

    Returns:
        dict: input_id : chembl_id dict
    """
    if not linker_fn:
        linker_fn = CHEMBL_DIR / 'unichem_linker.json'
    if not isinstance(x, list):
        x = [x]
        x = flatten(x)
    if id_from == 'cid':
        chembl_dict = load_chembl_pc_linker(linker_fn)
        ret_dict = {}
        for c in x:
            try:
                ret_dict[c] = chembl_dict[c]
            except KeyError as ke:
                ret_dict[c] = None
        return ret_dict
    elif id_from == 'uniprot':
        targ = new_client.target
        res = targ.filter(target_components__accession__in=x)
        ret_dict = {i['target_components'][0]['accession'] : i['target_chembl_id'] for i in res}
        return ret_dict
    else:
        raise ValueError(f'Unsupported id_from {id_from}')


def get_uniprot_acc(x, id_from='entry_name'):
    """Convert input  to uniprot accession

    Currently only supports Uniprot entry name conversion

    Args:
        x (str or list): Values to be converted
        id_from (str, optional): Value type to be converted. Defaults to 'entry_name'.

    Raises:
        ValueError: raised if unsupported id_from provided

    Returns:
        dict: input : uniprot_accession dict
    """
    if not isinstance(x, list):
        x = [x]
    if id_from == 'entry_name':
        chunks = list_chunk(x, 50)
        resp_df = pd.DataFrame()
        for chunk in tqdm(chunks, bar_format='{n_fmt} / {total_fmt}', total=ceil(len(x) / 50)):
            if isinstance(chunk[0],list): # multiple uniprot ids per entry
                chunk = [' '.join(c) for c in chunk]
            params = {
            'from': 'ACC+ID',
            'to': 'ACC',
            'format': 'tab', # for some probably masochistic reason, uniprot doesn't support json return format
            'query': ' '.join(chunk)
            }
            data = parse.urlencode(params)
            data = data.encode('utf-8')
            req = request.Request(UNIPROT_URL, data)
            with request.urlopen(req) as f:
                response = f.read()
                resp_dat = response.decode('utf-8')
            resp_df = resp_df.append(pd.read_csv(StringIO(resp_dat), sep='\t'))
            sleep(0.05)
        resp_dict = dict(zip(map(str, resp_df['From']), map(str, resp_df['To'])))
        return resp_dict
    else:
        raise ValueError(f'Unsupported id_from {id_from}')


def get_cid(x, id_from='sid'):
    """Convert input to pubchem CID

    Currently only supports conversion from pubchem SID

    Args:
        x (str or list): input to be converted
        id_from (str, optional): input type to be converted. Defaults to 'sid'.

    Raises:
        ValueError: raised if unsupported id_from provided

    Returns:
        dict: input : pubchem_cid dict
    """
    if not isinstance(x, list):
        x = [x]
    if id_from == 'sid':
        cs_dict = {}
        for chunk in tqdm(list_chunk(x, 100), bar_format='{n_fmt} / {total_fmt}', total=ceil(len(x) / 100)):
            resp = request.urlopen(PUG_REST_URL.format(','.join(chunk)))
            data = json.loads(resp.read())
            for dct in data['InformationList']['Information']:
                if 'CID' in dct.keys():
                    cs_dict.update({str(dct['SID']) : str(dct['CID'][0])})
                else:
                    cs_dict.update({str(dct['SID']) : None})
            sleep(0.05)
        return(cs_dict)
    else:
        raise ValueError(f'Unsupported id_from {id_from}')


