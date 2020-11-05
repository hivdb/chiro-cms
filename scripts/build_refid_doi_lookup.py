#! /usr/bin/env python3
import sys
import csv
import json
import requests
from more_itertools import chunked

IDCONV_API = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/'
IDCONV_EMAIL = 'hivdbteam@stanford.edu'


def rows2lookup(rows):
    rows = list(rows)
    pmids = [row['PMID'].strip() for row in rows]
    pmids = [pmid for pmid in pmids if pmid and pmid.isdigit()]
    pmid2doi = {}
    for partial in chunked(pmids, 150):
        resp = requests.post(
            IDCONV_API,
            {'tool': 'hivdb-pmid2doi',
             'email': IDCONV_EMAIL,
             'ids': ','.join(partial),
             'format': 'json'})
        data = resp.json()
        for one in data['records']:
            if 'doi' in one:
                pmid2doi[one['pmid']] = one['doi']
    results = []
    for row in rows:
        refid = row['RefID'].strip()
        extid = row['PMID'].strip()
        if extid in pmid2doi:
            extid = pmid2doi[extid]
        if extid.startswith('10.'):
            results.append({
                'refId': refid,
                'doi': extid
            })
    return results


def main():
    if len(sys.argv) != 3:
        print('Usage: {} <INPUT_CSV> <OUTPUT_JSON>', file=sys.stderr)
        exit(1)
    in_csv, out_json = sys.argv[1:]
    results = []
    with open(in_csv, encoding='utf-8-sig') as fp:
        reader = csv.DictReader(fp)
        results = rows2lookup(reader)
    with open(out_json, 'w') as fp:
        json.dump(results, fp, indent=2)


if __name__ == '__main__':
    main()
