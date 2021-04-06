#! /usr/bin/env python3
import sys
import csv
import json
import requests
from more_itertools import chunked

ENTREZ_API = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
NCBI_APIKEY = 'c589b6589a876ae42089c059c49249722807'


def rows2lookup(rows):
    rows = list(rows)
    pmids = [
        row['DOI'].strip()
        if row['DOI'].strip()
        else row['PMID'].strip()
        for row in rows
    ]
    pmids = [pmid for pmid in pmids if pmid and pmid.isdigit()]
    pmid2doi = {}
    for partial in chunked(pmids, 150):
        resp = requests.post(
            ENTREZ_API,
            {'db': 'pubmed',
             'id': ','.join(partial),
             'format': 'json',
             'api_key': NCBI_APIKEY})
        data = resp.json()
        pmids = data['result']['uids']
        for pmid in pmids:
            one = data['result'][pmid]
            if 'error' in one:
                print('Error: {} ({})'.format(one['error'], pmid),
                      file=sys.stderr)
                exit(1)
            for extid in one['articleids']:
                if extid['idtype'] == 'doi':
                    pmid2doi[pmid] = extid['value']
                    break
    results = []
    for row in rows:
        refid = row['RefID'].strip()
        extid = row['DOI'].strip()
        if not extid:
            extid = row['PMID'].strip()
        if extid in pmid2doi:
            extid = pmid2doi[extid]
        if extid.startswith('10.'):
            results.append({
                'refId': refid,
                'doi': extid
            })
    results.sort(key=lambda row: row['refId'])
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
