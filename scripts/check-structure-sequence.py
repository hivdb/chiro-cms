import click
import csv
import re
from collections import defaultdict
from pathlib import Path


def _load_csv(file_path):
    records = []
    with open(file_path, encoding='utf-8-sig') as fd:
        for record in csv.DictReader(fd):
            records.append(record)

    return records


def load_csv(file_path):
    records = _load_csv(file_path)
    for row in records:
        for key in row.keys():
            value = row[key].replace('_x000D_\n', '\n')
            row[key] = value

    return records


def get_mab_from_mabs(mabs):
    mab_names = set()
    for item in mabs:
        mab_list = item['mAb (EC50)']
        mab_list = get_mab_list(mab_list)

        for name in mab_list:
            mab_names.add(name)

    mab_names = list(mab_names)
    mab_names.sort()

    return mab_names


def get_mab_list(mab_list):
    mab_list = re.split('(?:[;\n]{1})', mab_list)
    mab_list = [i.strip() for i in mab_list]

    result = []
    for name in mab_list:
        if '(' in name:
            name = name.split('(', 1)[0]
        name = name.strip()
        if name:
            result.append(name)

    return result


def get_mab2refIDs_from_mabs(mabs):
    mab2refID = defaultdict(list)
    for item in mabs:
        mab_list = item['mAb (EC50)']
        mab_list = get_mab_list(mab_list)

        refIDs = item['RefIDs']
        refIDs = re.split('(?:[;\n]{1})', refIDs)
        refIDs = [i.strip() for i in refIDs]
        refIDs = [i for i in refIDs if i]
        refIDs = ';'.join(refIDs)
        if not refIDs:
            refIDs = 'NA'

        for name in mab_list:
            mab2refID[name].append(refIDs)

    return mab2refID


def get_mab_from_structures(structures):
    mab_names = set()
    for item in structures:
        name = item['mAb Names']
        name = name.strip()
        if ';' in name:
            continue

        if name:
            mab_names.add(name)

    mab_names = list(mab_names)
    mab_names.sort()

    return mab_names


def get_mab2refIDs_from_structures(structures):
    mab2refID = defaultdict(list)
    for item in structures:
        name = item['mAb Names']
        name = name.strip()
        if ';' in name:
            continue

        refIDs = item['RefIDs']
        if name:
            mab2refID[name] = refIDs

    return mab2refID


def validate_structures(structures):

    for item in structures:
        if not item['RefIDs']:
            print('Structure Error, no RefIDs', item)

    for item in structures:
        refIDs = item['RefIDs']
        if refIDs == 'NA':
            continue
        if '20' not in refIDs:
            print('RefID dont have publish year', refIDs)


def get_mab_has_structures(structures):
    mab_names = set()
    for item in structures:
        name = item['mAb Names']
        name = name.strip()
        if ';' in name:
            continue

        epitope = item['Epitope(<=4.5Å)']
        if not epitope:
            continue
        if name:
            mab_names.add(name)

    mab_names = list(mab_names)
    mab_names.sort()

    return mab_names


def get_mab_from_sequences(sequences):
    mab_names = set()
    for item in sequences:
        name = item['Mab']
        name = name.strip()
        if name:
            mab_names.add(name)

    mab_names = list(mab_names)
    mab_names.sort()

    return mab_names


def get_mab2refIDs_from_sequences(sequences):
    mab2refID = defaultdict(list)
    for item in sequences:
        name = item['Mab']
        name = name.strip()
        refIDs = item['Author']
        if name:
            mab2refID[name] = refIDs

    return mab2refID


def validate_sequences(sequences):
    for item in sequences:
        refIDs = item['Author']
        if refIDs == 'NA':
            continue
        if '20' not in refIDs:
            print('RefID dont have publish year', refIDs)


def show_mab_diff(list1, list1_name, list2, list2_name):
    diff = set(list1) - set(list2)
    print('In {} not in {}:'.format(list1_name, list2_name))
    for i in sorted(diff):
        print(i)
    print('=' * 20)


def show_error_refID(base, mapper):
    for k, v in mapper.items():
        base_v = base[k]
        if not base_v:
            continue

        if v not in base_v:
            print(k, v, base_v)

    print('=' * 20)


@click.command()
@click.argument('folder_path', type=click.Path(exists=True))
def work(folder_path):
    folder = Path(folder_path).absolute()

    mabs = load_csv(folder / 'MAbs.csv')
    structures = load_csv(folder / 'Structures.csv')
    sequences = load_csv(folder / 'Sequences.csv')

    print('Verify structures')
    validate_structures(structures)
    print('Verify sequences')
    validate_sequences(sequences)

    mab_names1 = get_mab_from_mabs(mabs)
    mab_names2 = get_mab_from_structures(structures)
    mab_names3 = get_mab_from_sequences(sequences)
    show_mab_diff(mab_names1, 'Mabs.csv', mab_names2, 'Structures.csv')
    print('Should be 0.')
    show_mab_diff(mab_names1, 'Mabs.csv', mab_names3, 'Sequences.csv')
    print('Should be 0.')

    mab2RefID_base = get_mab2refIDs_from_mabs(mabs)
    mab2RefID_struct = get_mab2refIDs_from_structures(structures)
    mab2RefID_seq = get_mab2refIDs_from_sequences(sequences)

    print('Structure RefIDs')
    show_error_refID(mab2RefID_base, mab2RefID_struct)
    print('Should be 0.')
    print('Sequence RefIDs')
    show_error_refID(mab2RefID_base, mab2RefID_seq)
    print('Should be 0.')

    mab_names2 = get_mab_has_structures(structures)
    show_mab_diff(mab_names2, 'Structures.csv', mab_names3, 'Sequences.csv')


if __name__ == '__main__':
    work()
