#! /usr/bin/env python

import csv
from pathlib import Path

import click
import ruamel.yaml


yaml = ruamel.yaml.YAML()


def inline_list(*lst):
    ret = ruamel.yaml.comments.CommentedSeq(lst)
    ret.fa.set_flow_style()
    return ret


def build_epitope_yaml_lookup(annot_dir):
    annot_dir.mkdir(parents=True, exist_ok=True)
    lookup = {}
    for yamlfile in annot_dir.iterdir():
        if yamlfile.suffix.lower() in ('.yml', 'yaml'):
            with yamlfile.open() as fp:
                data = yaml.load(fp)
                if isinstance(data, dict):
                    data = [data]
                if not data:
                    raise RuntimeError('File {} is empty'.format(yamlfile))
                for one in data:
                    name = one['name']
                    lookup[name] = yamlfile
    return lookup


def load_epitope_csv(fp):
    epitope_pairs = []
    for row in csv.DictReader(fp):
        pair = {}
        for key, val in row.items():
            lower_key = key.lower()
            if lower_key in ('mab', 'mab names'):
                pair['mab'] = val.strip()
            elif lower_key.startswith('epitope'):
                pair['epitopes'] = [
                    int(e) for e in val.strip().split('+') if e
                ]
            elif lower_key == 'pdb' and val:
                pair['pdb'] = val
        if pair['epitopes']:
            epitope_pairs.append(pair)
    return epitope_pairs


def save_epitopes(mab, epitopes, pdb, dest_yaml):
    if dest_yaml.is_file():
        with dest_yaml.open() as fp:
            payload = yaml.load(fp)
    else:
        payload = {
            'name': mab,
            'level': 'position',
            'hideCitations': True,
            'colorRules': [],
            'citations': [{
                'doi': '!!!ADD_BEFORE_COMMIT!!!',
                'author': '!!!ADD_BEFORE_COMMIT!!!',
                'year': '!!!ADD_BEFORE_COMMIT!!!',
                'section': 'PDB: {}'.format(pdb),
                'positions': []
            }],
            'positions': []
        }
    if isinstance(payload, list):
        target = [p for p in payload if p['name'] == mab][0]
    else:
        target = payload
    pos_target = [p for p in target['positions'] if p.get('subgroup') == mab]
    if pos_target:
        pos_target = pos_target[0]
    else:
        pos_target = {
            'subgroup': mab,
            'positions': []
        }
        target['positions'].append(pos_target)
    if (*pos_target['positions'],) == (*epitopes,):
        # no change is necessary
        return
    pos_target['positions'] = inline_list(*epitopes)
    target['citations'][0]['positions'] = inline_list(*epitopes)
    if dest_yaml.is_file():
        print('Updated {}, manual review is required'.format(dest_yaml))
    else:
        print('Added {}, manual review is required'.format(dest_yaml))

    with dest_yaml.open('w') as fp:
        yaml.dump(payload, fp)


@click.command()
@click.argument('input_csv',
                type=click.File('r', encoding='utf-8-sig'))
@click.argument('annot_directory',
                type=click.Path(dir_okay=True, file_okay=False))
def import_epitopes(input_csv, annot_directory):
    epitope_pairs = load_epitope_csv(input_csv)
    annot_directory = Path(annot_directory)
    yaml_lookup = build_epitope_yaml_lookup(annot_directory)
    for pair in epitope_pairs:
        mab = pair['mab']
        if mab in yaml_lookup:
            dest_yaml = yaml_lookup[mab]
        else:
            dest_yaml = annot_directory / '{}.yml'.format(mab.lower())
        save_epitopes(dest_yaml=dest_yaml, **pair)


if __name__ == '__main__':
    import_epitopes()
