#! /usr/bin/env python

import os
import re
import csv
import sys
import textwrap
import requests
import ruamel.yaml
from ruamel.yaml.scalarstring import PreservedScalarString

from functools import lru_cache


yaml = ruamel.yaml.YAML()
yaml.indent(sequence=4, offset=2)


def wrapped(s, width=70):
    if '\n' in s or '\r' in s:
        return PreservedScalarString(s)
    return PreservedScalarString('\n'.join(textwrap.wrap(s, width=width)))


def norm_space(text):
    if text is None:
        return None
    return re.sub(
        r'[ \t]+[\r\n]+', '\n',
        re.sub(r'[ \t]+', ' ', text)
    )


def get_first(row, *keys, required=True, type_coerce=str):
    if not keys:
        raise TypeError(
            'get_first() argument "keys" can not be empty'
        )
    for key in keys:
        val = row.get(key)
        if val is not None:
            val = norm_space(val.strip())
            if val == 'NA':
                return None
            return type_coerce(val)
    if required and val is None:
        raise KeyError(
            'CSV file don\'t have "{}" column'
            .format(keys[0])
        )


def smart_split2(text, seps=';\r\n'):
    pattern = re.compile(r'[{}]+'.format(re.escape(seps)))
    if not text:
        return []
    parts = [r.strip() for r in pattern.split(text)]
    return [r for r in parts if r]


DOI_PATTERN = re.compile(r'(10\.\d{4,}(?:\.\d+)*/[^)\s]+)')


@lru_cache
def read_doi(doi):
    resp = requests.get(
        'https://api.crossref.org/works/{}'.format(doi)
    )
    resp.raise_for_status()
    data = resp.json()
    return data['message']


def extract_first_author(result):
    for author in result.get('author', []):
        if 'family' not in author:
            continue
        return {
            'surname': author['family'],
            'givenNames': author.get('given')
        }


def extract_year(result):
    years = []
    for key in ('published-print', 'created', 'issued'):
        if key in result and 'date-parts' in result[key]:
            year = result[key]['date-parts'][0][0]
            if year:
                years.append(year)
    return min(years) if years else None


def extract_journal_title(result):
    for jtitle in result.get('container-title', []):
        return jtitle
    if 'institution' in result:
        return result['institution']['name']


def extract_short_journal_title(result):
    for jtitle in result.get('short-container-title', []):
        return jtitle


def extract_references(row):
    refs = get_first(row, 'references')
    if not refs:
        return []
    refobjs = []
    for doi in DOI_PATTERN.findall(refs):
        message = read_doi(doi)
        journal = extract_journal_title(message)
        journal_short = extract_short_journal_title(message)

        refobj = {
            'doi': doi,
            'firstAuthor': extract_first_author(message),
            'year': extract_year(message),
            'journal': journal
        }
        if journal_short and journal != journal_short:
            refobj['journalShort'] = journal_short
        refobjs.append(refobj)
    if not refobjs:
        raise ValueError("No DOI was found in row {}"
                         .format(get_first(row, 'mab names')))
    return refobjs


def extract_source(row):
    return smart_split2(get_first(row, 'source', 'sources'))


EC50_PATTERN = re.compile(r"""
  \s*
  (>|<|=|<=|>=|~|<<)?\s*           # cmp
  (-?\d[\d,]*\.?(?:\d+)?)\s*    # num
  ([\xb5\u03bcu]M|IU/ml|U/ml|ng/ml)?   # unit
  .*
""", re.X | re.I)


unit_variants = {
    'um': '\xb5M',
    '\xb5m': '\xb5M',
    '\u03bcm': '\xb5M',  # greek letter \u03bc to micro sign \xb5
    'iu/ml': 'IU/ml',
    'u/ml': 'U/ml',
    'ng/ml': 'ng/ml',
}


def contains_word(text, word):
    match = re.search(r'\b{}\b'.format(re.escape(word)), text)
    return bool(match)


def extract_ec50(row, ab_names):
    ec50s = smart_split2(get_first(row, 'ec50'))
    ec50lookup = {}
    for ec50desc in ec50s:
        for ab_name in ab_names:
            if not contains_word(ec50desc, ab_name):
                continue
            ec50obj = {}
            for part in ec50desc.split(ab_name):
                match = EC50_PATTERN.search(part)
                if match:
                    result = []
                    unit = match.group(3) or 'ng/ml'
                    unit = unit_variants[unit.lower()]
                    ec50cmp = match.group(1) or '='
                    number = float(match.group(2).replace(',', ''))
                    if ec50cmp != '=':
                        result.append(ec50cmp)
                    result.append('{:2f}'.format(number)
                                  .rstrip('0').rstrip('.'))
                    if unit != 'ng/ml':
                        result.append(unit)
                    if len(result) == 1:
                        result = result[0]
                        if '.' in result:
                            result = float(result)
                        else:
                            result = int(result)
                    else:
                        result = ''.join(result)
                    ec50obj['ec50'] = result
                if 'PV' in part:
                    ec50obj['ec50Note'] = 'PV'
                elif 'IC100' in part:
                    ec50obj['ec50Note'] = 'IC100'
                elif 'IC90' in part:
                    ec50obj['ec50Note'] = 'IC90'
            ec50lookup[ab_name] = ec50obj
            break
    return ec50lookup


def extract_pdb(row, ab_names):
    pdbs = smart_split2(get_first(row, 'pdb structures', 'pdb'))
    pdblookup = {}
    for pdbdesc in pdbs:
        for ab_name in ab_names:
            if not contains_word(pdbdesc, ab_name):
                continue
            for part in pdbdesc.split(ab_name):
                part = part.strip(': \t')
                if part:
                    pdblookup[ab_name] = part
                    break
            break
    return pdblookup


def extract_antibodies(row, alnlookup):
    ab_names = smart_split2(get_first(
        row, 'mab names', 'ab names', 'antibody names'))
    ec50lookup = extract_ec50(row, ab_names)
    pdblookup = extract_pdb(row, ab_names)
    ab_objs = []
    for ab_name in ab_names:
        ab_obj = {
            'name': ab_name
        }
        ec50obj = ec50lookup.get(ab_name)
        if ec50obj:
            ab_obj.update(ec50obj)
        pdbdesc = pdblookup.get(ab_name)
        if pdbdesc:
            ab_obj['pdb'] = pdbdesc
        alndata = alnlookup.get(ab_name)
        if alndata:
            ab_obj.update(alndata)
        ab_objs.append(ab_obj)
    return sorted(ab_objs, key=lambda obj: obj['name'])


def load_data(filepath):
    with open(filepath, encoding='utf-8-sig') as fp:
        rows = []
        for row in csv.DictReader(fp):
            lowerkey_row = {}
            for key in row:
                lowerkey_row[key.lower()] = row[key]
            rows.append(lowerkey_row)
        return rows


def load_aln_data(filepath):
    with open(filepath, encoding='utf-8-sig') as fp:
        lookup = {}
        for row in csv.DictReader(fp):
            ab_name = row['Mab']
            cdrh3len = row['CDRH3 Length']
            cdrl3len = row['CDRL3 Length']
            cdrh3len = int(cdrh3len) if cdrh3len else None
            cdrl3len = int(cdrl3len) if cdrl3len else None
            lookup[ab_name] = {
                'species': row['Species'],
                'IGHV': row['IGHV'],
                'PcntMutH': row['%Mut(H)'],
                'CDRH3Len': cdrh3len,
                'IGHJ': row['IGHJ'],
                'IGLV': row['IGLV'],
                'PcntMutL': row['%Mut(L)'],
                'CDRL3Len': cdrl3len,
                'IGLJ': row['IGLJ']
            }
        return lookup


def yaml_filename(groupobj):
    mabnames = []
    for mab in groupobj['antibodies']:
        mabnames.append(mab['name'])
    return '_'.join(mabnames) + '.yml'


def main():
    if len(sys.argv) != 4:
        print('Usage: {} <MAbCSV> <MAbAlignmentCSV> <OUTDIR>'
              .format(sys.argv[0]), file=sys.stderr)
        exit(1)
    mabcsv_path, alncsv_path, outdir = sys.argv[1:]
    os.makedirs(outdir, exist_ok=True)

    mabdata = load_data(mabcsv_path)
    alnlookup = load_aln_data(alncsv_path)
    for row in mabdata:
        groupobj = {
            'references': extract_references(row),
            'source': extract_source(row),
            'antibodies': extract_antibodies(row, alnlookup),
            'config': wrapped(
                get_first(row, 'configuration', 'config')
            )
        }
        with open(
            os.path.join(
              outdir, yaml_filename(groupobj)
            ), 'w'
        ) as fp:
            yaml.dump(groupobj, fp)


if __name__ == '__main__':
    main()
