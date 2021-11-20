import os
import sys
import requests
import re
import time

BASE_URL = 'https://api.fda.gov/drug/ndc.json?'
CLASS_TERM = '&search=openfda.pharm_class_epc:""'
GENERIC_NAME = '&search=generic_name:'
BRAND_NAME = '&search=brand_name:' 
API_BASE = "api_key="
LIMIT_1 = '&limit=1'

epc_regex = re.compile("EPC")

def get_pharm_classes(terms, all_classes = True, api_key=None, dedup=True):
    terms = [terms] if isinstance(terms, str) else terms
    classes = {}
    for t in terms:
        if api_key:
            gen_resp = requests.get(url = BASE_URL+API_BASE+api_key+GENERIC_NAME+t+LIMIT_1)
        else:
            gen_resp = requests.get(url = BASE_URL+GENERIC_NAME+t+LIMIT_1)
        gen_dat = gen_resp.json()
        try:
            pharm_classes = gen_dat['results'][0]['pharm_class']
        except KeyError:
            pharm_classes = "N/A"
        finally:
            if dedup:
                pharm_classes = list(set(pharm_classes))
            if all_classes:
                classes[t] = pharm_classes
            else:
                classes[t] = [m for m in pharm_classes if epc_regex.search(m)]
            time.sleep(0.1)
    return classes


if __name__ == "__main__":
    res = get_pharm_classes(["atorvastatin calcium", "aspirin"], all_classes = False)
    print(res)


