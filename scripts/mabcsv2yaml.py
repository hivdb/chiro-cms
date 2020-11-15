#! /usr/bin/env python

import re
import csv
import json
from functools import lru_cache

import click
import requests
import ruamel.yaml
from ruamel.yaml.scalarstring import PreservedScalarString


yaml = ruamel.yaml.YAML()


def norm_space(text):
    if text is None:
        return None
    return re.sub(
        r'[ \t]+[\r\n]+', '\n',
        re.sub(r'[ \t]+', ' ', text)
    )


def get_first(row, *keys, required=True, type_coerce=str, null_value='NA'):
    if not keys:
        raise TypeError(
            'get_first() argument "keys" can not be empty'
        )
    for key in keys:
        val = row.get(key)
        if val is not None:
            val = norm_space(val.strip())
            if val == null_value:
                return None
            return type_coerce(val)
    if required and val is None:
        raise KeyError(
            'CSV file don\'t have "{}" column'
            .format(keys[0])
        )


def smart_split2(text, seps=';\r\n'):
    if not text:
        return []
    pattern = re.compile(r'[{}]+'.format(re.escape(seps)))
    parts = [r.strip() for r in pattern.split(text)]
    return [r for r in parts if r]


REFID_PATTERN = re.compile(r'\b[\w-]+\d{2}[a-z]?\b')


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
    if 'publisher' in result:
        return result['publisher']


def extract_short_journal_title(result):
    for jtitle in result.get('short-container-title', []):
        return jtitle


def extract_references(row, refid_lookup):
    refids = get_first(row, 'refids')
    free_text_refs = get_first(row, 'freetextrefs', required=False)
    if not refids and not free_text_refs:
        return []
    refobjs = []
    for refid in REFID_PATTERN.findall(refids):
        doi = refid_lookup.get(refid)
        if not doi:
            raise click.ClickException('Unable to find RefID {}'.format(refid))
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
    for free_text in smart_split2(free_text_refs):
        refobjs.append({'freeText': free_text})
    return refobjs


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
    match = re.search(r'(?<![\w\'"`+-]){}(?![\w\'"`+-])'
                      .format(re.escape(word)), text)
    return bool(match)


def parse_ec50(ec50desc):
    ec50obj = {}
    match = EC50_PATTERN.search(ec50desc)
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
    if 'PV' in ec50desc:
        ec50obj['ec50Note'] = 'PV'
    elif 'IC100' in ec50desc:
        ec50obj['ec50Note'] = 'IC100'
    elif 'IC90' in ec50desc:
        ec50obj['ec50Note'] = 'IC90'
    return ec50obj


def extract_antibodies(row):
    antibodies = get_first(row, 'mab (ec50)', required=False)
    if not antibodies:
        return {}
    result = []
    for antibody in smart_split2(antibodies):
        antibody = antibody.split(' ', 1)
        if len(antibody) == 1:
            result.append({
                'name': antibody[0]
            })
        else:
            antibody, ec50desc = antibody
            ec50obj = parse_ec50(ec50desc)
            result.append({
                'name': antibody,
                **ec50obj
            })
    return result


def load_data(fp):
    rows = []
    for row in csv.DictReader(fp):
        lowerkey_row = {}
        for key in row:
            row[key] = row[key].replace('_x000D_\n', '\n')
            lowerkey_row[key.lower()] = row[key]
        rows.append(lowerkey_row)
    return rows


def extract_igxx(text, title, pcnt_title=None):
    text = text.strip()
    result = {title: None}
    if pcnt_title:
        result[pcnt_title] = None
    if not text:
        return result
    match = re.match(r'^([^*]+)(?:\*[^(]+(?:\((\d+(?:\.\d+)?%)\))?)?$',
                     text.strip())
    if not match:
        raise click.ClickException('Invalid IGXX string: {!r}'.format(text))
    groups = match.groups()
    result[title] = re.sub(r'^IG', '', groups[0])
    result[title] = result[title].replace(' (', '\n(')
    if pcnt_title and len(groups) == 2:
        result[pcnt_title] = groups[1]
    return result


def extract_cdrlen(text):
    text = text.strip()
    match = re.search(r'\((\d+)\)$', text)
    if match:
        return int(match.group(1))


def load_refid_lookup(fp):
    pair = json.load(fp)
    lookup = {}
    for row in pair:
        lookup[row['refId']] = row['doi']
    return lookup


def yaml_filename(groupobj):
    mabnames = []
    for mab in groupobj['antibodies']:
        mabnames.append(mab['name'])
    return '_'.join(mabnames) + '.yml'


@click.command()
@click.argument('input_mabcsv',
                type=click.File('r', encoding='utf-8-sig'))
@click.option(
    '--refid-lookup', type=click.File('r', encoding='utf-8-sig'),
    help='RefID2DOI lookup JSON file')
@click.argument('output_yaml', type=click.File('w'))
def mabcsv2yaml(input_mabcsv, refid_lookup, output_yaml):

    mabdata = load_data(input_mabcsv)
    refid_lookup = load_refid_lookup(refid_lookup)
    groups = []
    for row in mabdata:
        groupobj = {
            'references': extract_references(row, refid_lookup),
            'antibodies': extract_antibodies(row),
        }
        source = get_first(row, 'source', 'sources')
        if source:
            groupobj['source'] = PreservedScalarString(source)

        desc = get_first(row, 'description')
        if desc:
            groupobj['description'] = PreservedScalarString(desc)

        groupobj['dataAvailability'] = {}
        if get_first(row, 'seq'):
            groupobj['dataAvailability']['sequence'] = True
        if get_first(row, 'pdb'):
            groupobj['dataAvailability']['structure'] = True
        if get_first(row, 'animal'):
            groupobj['dataAvailability']['animal'] = True
        if get_first(row, 'trial'):
            groupobj['dataAvailability']['trial'] = True

        groups.append(groupobj)
    yaml.dump(groups, output_yaml)


if __name__ == '__main__':
    mabcsv2yaml()
