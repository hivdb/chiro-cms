#! /usr/bin/env python

import re
import sys
import ruamel.yaml
from pathlib import Path
from itertools import chain
from collections import defaultdict


yaml = ruamel.yaml.YAML()


def parse_mut(mut):
    match = re.match(
        r'^[A-Z]?(\d+)([A-Zid*]+|[A-Z*]_[A-Z*]+)$', mut)
    if not match:
        raise ValueError('Invalid mutation string: {!r}'.format(mut))
    pos, aas = match.groups()
    pos = int(pos)
    if '_' in aas:
        aas = {aas}
    else:
        aas = set(aas)
    return pos, aas


def parse_muts(muts):
    posaas = defaultdict(set)
    for mut in muts:
        pos, aas = parse_mut(mut)
        posaas[pos] |= aas
    return posaas


def compose_mut(pos, aas):
    muts = []
    remains = []
    for aa in aas:
        if '_' in aa:
            muts.append('{}{}'.format(pos, aa))
        remains.append(aa)
    muts.append('{}{}'.format(pos, ''.join(sorted(remains))))
    return sorted(muts)


def compose_muts(posaas):
    muts = []
    for pos, aas in sorted(posaas.items()):
        muts.extend(compose_mut(pos, aas))
    return muts


def convert(mab, byref_data, bymab_data):
    if bymab_data is None:
        bymab_data = {
            'name': 'resistance {}'.format(mab),
            'label': mab,
            'tableColLabel': '{} DRM'.format(mab),
            'level': 'aminoAcid',
            'hideCitations': False,
            'colorRules': [],
            'citations': [],
            'aminoAcids': []
        }
    posaas = parse_muts(bymab_data['aminoAcids'])
    cite_posaas = defaultdict(set)
    for mut, mabs in byref_data['aminoAcidAttrs'][0]['pairs'].items():
        pos, aas = parse_mut(mut)
        if mab in mabs:
            posaas[pos] |= aas
            cite_posaas[pos] |= aas
    bymab_data['citations'].append({
        'doi': byref_data['citations'][0]['doi'],
        'author': byref_data['citations'][0]['author'],
        'year': byref_data['citations'][0]['year'],
        'section': byref_data['citations'][0]['section'],
        'positions': sorted(cite_posaas.keys()),
        'aminoAcids': compose_muts(cite_posaas)
    })
    bymab_data['aminoAcids'] = compose_muts(posaas)
    return bymab_data


def main():
    destdir = Path(sys.argv[1])
    all_bymab_data = {}
    for ymlfile in destdir.glob('escape-mutants-*.yml'):
        with ymlfile.open() as fp:
            byref_data = yaml.load(fp)
            mabs = set(chain(*list(
                byref_data['aminoAcidAttrs'][0]['pairs'].values()
            )))
            for mab in mabs:
                bymab_data = all_bymab_data.get(mab)
                bymab_data = convert(mab, byref_data, bymab_data)
                all_bymab_data[mab] = bymab_data
        ymlfile.unlink()
    for mab, bymab_data in all_bymab_data.items():
        destfile = destdir / 'resistance-{}.yml'.format(
            mab.lower().replace('+', 'p')
        )
        with open(destfile, 'w') as fp:
            yaml.dump(bymab_data, fp)


if __name__ == '__main__':
    main()
