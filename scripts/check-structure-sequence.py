import click
import csv
import re
from collections import defaultdict
from pathlib import Path

ABDAB_FILE = 'CoV-AbDab_261120.csv'
RENAME_MAB_ABDAB2MAb = {
    'MR17_K99Y': 'MR17-K99Y',
    'S2-M11': 'S2M11',
    'S2-E12': 'S2E12',
    'Fab2-4': '2-4',
    'BD23': 'BD-23',
    'Nb20': 'nb20',
    'Nb#6': 'Nb6',
}


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


def strip_then_remove(a_list):
    """caller must make sure all items are string"""
    a_list = [i.strip() for i in a_list]
    a_list = [i for i in a_list if i]
    return a_list


def load_abdab(folder):
    file_path = folder / ABDAB_FILE
    return load_csv(file_path)


def get_mab2PDB_from_abdab(abdab):
    mab2PDB = defaultdict(list)
    for item in abdab:
        if not is_sars_cov2_mab(item):
            continue

        name = item['Name']
        name = name.strip()
        pdb = item['Structures']
        pdb = pdb.strip()
        if pdb in ['ND', 'MD']:
            continue

        if pdb.startswith('Expected'):
            pdb = process_abdab_expected(pdb)
            if pdb:
                mab2PDB[name].extend(pdb)
        elif pdb.startswith('http'):
            if 'Expected' in pdb:
                expected = pdb[pdb.index('Expected'):]
                pdb = pdb[:pdb.index('Expected')]
            else:
                expected = None
            pdb = re.split(r'(?:[;,]{1})', pdb)
            pdb = [re.sub(
                r'https?://www.rcsb.org/structure/', '', i) for i in pdb]
            pdb = strip_then_remove(pdb)
            if not pdb and not expected:
                continue
            if expected:
                expected = process_abdab_expected(expected)
                pdb.extend(expected)
            mab2PDB[name].extend(pdb)
        elif pdb.startswith('ND'):
            print('Ignore:', name, pdb)
        elif pdb.startswith('See'):
            print('Ignore:', name, pdb)
        else:
            mab2PDB[name].append(pdb)

    for k, v in mab2PDB.items():
        for i in v:
            if i.startswith('!'):
                i = i[1:]
            if len(i) != 4:
                print(k, v)

    return mab2PDB


def get_mab2Seq_from_abdab(abdab):
    mab2Seq = defaultdict(dict)

    for item in abdab:
        if not is_sars_cov2_mab(item):
            continue

        name = item['Name']
        name = name.strip()

        if name == 'Nb6':
            # TODO, RefID + mab
            continue

        if name in RENAME_MAB_ABDAB2MAb:
            name = RENAME_MAB_ABDAB2MAb[name]

        vh = item['VH or VHH'].strip()
        if vh in ['N/A', 'ND', 'Expected']:
            vh = ''
        vl = item['VL'].strip()
        if vl in ['N/A', 'ND', 'Expected']:
            vl = ''
        cdrh3 = item['CDRH3'].strip()
        if cdrh3 in ['N/A', 'ND', 'Expected']:
            cdrh3 = ''
        cdrl3 = item['CDRL3'].strip()
        if cdrl3 in ['N/A', 'ND', 'Expected']:
            cdrl3 = ''

        if name in mab2Seq:
            print('Duplicated name:', name)
        mab2Seq[name] = {
            'VH': vh,
            'VL': vl,
            'CDRH3': cdrh3,
            'CDRL3': cdrl3,
        }

    return mab2Seq


def process_abdab_expected(pdb):
    pdb = pdb.replace('Expected', '')
    pdb = pdb.strip()
    pdb = pdb[1:-1]
    pdb = re.split(r'(?:[;,]{1})', pdb)
    pdb = strip_then_remove(pdb)
    pdb = ['!{}'.format(i) for i in pdb]
    return pdb


def all_expected_pdb(pdb_list):
    result = [i.startswith('!') for i in pdb_list]
    if not result:
        return False
    return all(result)


def trans_expected_pdb(pdb_list):
    result = []
    for i in pdb_list:
        if i.startswith('!'):
            result.append(i[1:].upper())
        else:
            result.append(i.upper())
    return result


def is_sars_cov2_mab(record):
    binds_to = record['Binds to']
    neutralising = record['Neutralising Vs']
    if 'SARS-CoV2' in binds_to:
        return True
    if 'SARS-CoV2' in neutralising:
        return True
    return False


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
    mab_list = re.split(r'(?:[;\n]{1})', mab_list)
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
        refIDs = re.split(r'(?:[;\n]{1})', refIDs)
        refIDs = strip_then_remove(refIDs)
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

    # TODO: check dup mAbs
    # TODO: check use same PDB


def get_mab_has_epitope(structures):
    """ PDB file is released"""
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


def get_mab_has_pdbs(structures):
    mab_names = set()
    for item in structures:
        name = item['mAb Names']
        name = name.strip()
        if ';' in name:
            continue

        pdb = item['PDB']
        if not pdb:
            continue
        if name:
            mab_names.add(name)

    mab_names = list(mab_names)
    mab_names.sort()

    return mab_names


def get_mab2existPDBs(structures):
    mab2existPDBs = defaultdict(list)
    for item in structures:
        name = item['mAb Names']
        name = name.strip()
        if not name:
            continue

        if ';' in name:
            continue

        pdb = item['PDB']
        if '(' in pdb:
            pdb = pdb[:pdb.index('(')]
        pdb = [pdb]
        other_pdbs = item['other PDBs'].split(';')
        other_pdbs = strip_then_remove(other_pdbs)
        pdb.extend(other_pdbs)
        if not pdb:
            continue
        mab2existPDBs[name].extend(pdb)

    return mab2existPDBs


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

    # TODO: duplicated


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


def show_missing_structure(structures, mab2PDB):
    mab2existPDBs = get_mab2existPDBs(structures)
    mab_with_epitope = get_mab_has_epitope(structures)

    miss_pdb = set()
    for item, check_pdbs in mab2PDB.items():
        check_pdbs = trans_expected_pdb(check_pdbs)
        if item in RENAME_MAB_ABDAB2MAb:
            item = RENAME_MAB_ABDAB2MAb[item]
        if item not in mab2existPDBs:
            miss_pdb.add(item)
        else:
            exist_pdbs = mab2existPDBs[item]
            if set(exist_pdbs) - set(check_pdbs):
                print('Too many {}:'.format(item),
                      exist_pdbs, ';'.join(check_pdbs))
            if set(check_pdbs) - set(exist_pdbs):
                print('Too little {}:'.format(item),
                      exist_pdbs, ';'.join(check_pdbs))

    print('All missing PDB:', miss_pdb)

    miss_epitope = set()
    for name, pdbs in mab2PDB.items():
        if all_expected_pdb(pdbs):
            continue
        if name not in mab_with_epitope:
            miss_epitope.add(name)
    print('PDB without epitope process:', miss_epitope)


def show_error_sequences(sequences, mab2Seq):

    for item in sequences:
        name = item['Mab']
        if name in mab2Seq:
            base_line = mab2Seq[name]
        else:
            continue

        vh = item['VH or VHH']
        vl = item['VL']
        cdrh3 = item['CDRH3']
        cdrl3 = item['CDRL3']

        if vh != base_line['VH']:
            print('Seq Error VH:', name, base_line['VH'])
        if vl != base_line['VL']:
            print('Seq Error VL:', name, base_line['VL'])
        if base_line['CDRH3'] not in cdrh3:
            print('Seq Error CDRH3:', name, base_line['CDRH3'])
        if base_line['CDRL3'] not in cdrl3:
            print('Seq Error CDRL3:', name, base_line['CDRL3'])


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

    mab_names2 = get_mab_has_epitope(structures)
    show_mab_diff(mab_names2, 'Structures.csv', mab_names3, 'Sequences.csv')

    abdab = load_abdab(folder)
    mab2PDB = get_mab2PDB_from_abdab(abdab)
    show_missing_structure(structures, mab2PDB)

    mab2Seq = get_mab2Seq_from_abdab(abdab)
    show_error_sequences(sequences, mab2Seq)


if __name__ == '__main__':
    work()
