import yaml
from pathlib import Path
import re
from collections import defaultdict


post_s_gene = [
    'ORF3a',
    'E',
    'M',
    'ORF6',
    'ORF7a',
    'ORF7b',
    'ORF8',
    'N',
    'ORF10'
]


def load_yaml(file_path):
    return yaml.load(open(file_path), Loader=yaml.Loader)


def dump_yaml(file_path, content):
    with open(file_path, 'w') as fd:
        yaml.dump(content, fd, default_flow_style=False)


def work(src_file_path, dst_file_path, var_group):
    var_info = load_yaml(Path(src_file_path))
    config = load_yaml(dst_file_path)

    s_muts = defaultdict(int)
    non_s_muts = defaultdict(int)

    for i in var_info:
        for m in i.keys():
            if 'name' in m:
                continue
            if ':' in m:
                non_s_muts[m] += 1
            else:
                s_muts[m] += 1

    rbd, ntd, ctd = process_s_muts(s_muts)
    pre_s, post_s = process_non_s_muts(non_s_muts)

    subtables = {
        f'table_{var_group}_rbd': rbd,
        f'table_{var_group}_ntd': ntd,
        f'table_{var_group}_ctd': ctd,
        f'table_{var_group}_pre_s': pre_s,
        f'table_{var_group}_post_s': post_s,
    }

    for s, d in subtables.items():
        if s not in config['tables']:
            config['tables'][s] = {
                'columnDefs': [],
                'data': {
                    '_resource': 'mutation-viewer/compare.yml'
                },
                'tableScrollStyle': {
                    'maxHeight': 'none'
                }
            }

        config['tables'][s]['columnDefs'] = d

    dump_yaml(dst_file_path, config)


def process_non_s_muts(non_s_muts):
    non_s_muts = [m for m, c in non_s_muts.items() if c >= 2]

    pre_s_muts = sorted([
        m
        for m in non_s_muts
        if m.split(':', 1)[0].upper() not in post_s_gene
    ])
    post_s_muts = sorted([
        m
        for m in non_s_muts
        if m.split(':', 1)[0].upper() in post_s_gene
    ])

    pre_s_list = [{
        'name': 'name',
        'label': 'Variant'
    }]

    for i in pre_s_muts:
        gene, m = i.split(':')
        m = re.match(r'(\D+)(\d+)', m)
        ref, pos = m.groups()
        pre_s_list.append({
            'name': i,
            'label': f'{gene}:{ref} {pos}'
        })

    post_s_list = [{
        'name': 'name',
        'label': 'Variant'
    }]

    for i in post_s_muts:
        gene, m = i.split(':')
        m = re.match(r'(\D+)(\d+)', m)
        ref, pos = m.groups()
        post_s_list.append({
            'name': i,
            'label': f'{gene}:{ref} {pos}'
        })

    return pre_s_list, post_s_list


def process_s_muts(s_muts):

    s_muts = [m for m, c in s_muts.items() if c >= 2]

    rbd_muts = set()
    ntd_muts = set()
    ctd_muts = set()

    for i in s_muts:
        if i.upper() == 'D614':
            continue
        m = re.match(r'(\D+)(\d+)', i)
        ref, pos = m.groups()
        pos = int(pos)
        if pos < 306:
            ntd_muts.add((i, ref, pos))
        elif pos > 534:
            ctd_muts.add((i, ref, pos))
        else:
            rbd_muts.add((i, ref, pos))

    rbd_muts = sorted(rbd_muts, key=lambda x: x[-1])
    ntd_muts = sorted(ntd_muts, key=lambda x: x[-1])
    ctd_muts = sorted(ctd_muts, key=lambda x: x[-1])

    rbd_list = [{
            'name': 'name',
            'label': 'Variant'
        }] + [
        {
            'name': name,
            'label': f'{ref} {pos}'
        }
        for name, ref, pos in rbd_muts
    ]
    ntd_list = [{
            'name': 'name',
            'label': 'Variant'
        }] + [
        {
            'name': name,
            'label': f'{ref} {pos}'
        }
        for name, ref, pos in ntd_muts
    ]
    ctd_list = [{
            'name': 'name',
            'label': 'Variant'
        }] + [
        {
            'name': name,
            'label': f'{ref} {pos}'
        }
        for name, ref, pos in ctd_muts
    ]

    return rbd_list, ntd_list, ctd_list


if __name__ == '__main__':

    work(
        'resources/mutation-viewer/voc_mut_align.yml',
        'pages/variants.yml',
        'voc')

    work(
        'resources/mutation-viewer/omicron_mut_align.yml',
        'pages/variants.yml',
        'omicron')
