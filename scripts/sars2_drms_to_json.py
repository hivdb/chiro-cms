#! /usr/bin/env python3
import csv
import json

# This TSV file is generated by
# https://github.com/hivdb/sierra-sars2-paper/blob/main/resistance-mutations.rmd
with open('downloads/resistance-mutations/latest.tsv') as fp:
    rows = []
    for row in csv.DictReader(fp, delimiter='\t'):
        row['position'] = int(row['position'])
        row['invivo_selected'] = row['invivo_selected'] == 'TRUE'
        rows.append(row)

with open('downloads/resistance-mutations/latest.json', 'w') as fp:
    json.dump({
        'MAB': [row for row in rows if row['gene'] == 'S'],
        '_3CLPI': [row for row in rows if row['gene'] == '_3CLpro'],
        'RdRPI': [row for row in rows if row['gene'] == 'RdRP']
    }, fp, indent=2)

print("Create: downloads/resistance-mutations/latest.json")
