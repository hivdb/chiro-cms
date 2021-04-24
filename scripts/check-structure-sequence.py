import click
import csv
import re
from collections import defaultdict
from pathlib import Path

ABDAB_FILE = 'CoV-AbDab_210421.csv'
RENAME_MAB_ABDAB2MAb = {
    'MR17_K99Y': 'MR17-K99Y',
    'S2-M11': 'S2M11',
    'S2-E12': 'S2E12',
    'Fab2-7': '2-7',
    'Fab2-4': '2-4',
    'Fab2-51': '2-51',
    'Fab1-57': '1-57',
    'Fab5-24': '5-24',
    'Fab1-87': '1-87',
    'Fab4-18': '4-18',
    'BD23': 'BD-23',
    'Nb20': 'nb20',
    'Nb#6': 'Nb6',
    'Fab2-15': '2-15',
    'Regdanvimab': 'CT-P59',
    'DH1050-1': 'DH1050',
    'VHH_U': 'VHH U',
    'VHH_V': 'VHH V',
    'VHH_E': 'VHH E',
    'VHH_W': 'VHH W',
}

SWITCH_MAB = {
    'C1A-B12': 'C1A-B3',
    "C1A-B3": 'C1A-B12',
}

RENAME_MAB = {
    'LY-CoV016': 'CB6',
    'LY-CoV555': 'Bamlanivimab',
    'COV2-2050-LALA-PG': 'COV2-2050',
    'COV2-3025-LALA-PG': 'COV2-3025',
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


def dump_csv(file_path, records, headers=[]):
    if not records:
        return
    if not headers and records:
        headers = records[0].keys()

    with open(file_path, 'w', encoding='utf-8-sig') as fd:
        writer = csv.DictWriter(fd, fieldnames=headers)
        writer.writeheader()
        writer.writerows(records)


def strip_then_remove(a_list):
    """caller must make sure all items are string"""
    a_list = [i.strip() for i in a_list]
    a_list = [i for i in a_list if i]
    return a_list


def get_dup_key(records, key):
    tester = defaultdict(int)
    for i in records:
        value = i[key]
        if not value:
            continue
        tester[value] += 1

    dup_key = [i for i, j in tester.items() if j > 1]
    return dup_key


def get_invalid_null_records(records, keys):
    result = []
    for item in records:
        values = [item.get(key) for key in keys]
        value_exists = [i for i in values if i not in [None, '']]
        if not value_exists:
            continue
        if len(value_exists) != len(keys):
            result.append(item)

    return result


def get_depend_invalid_null_records(records, keys):
    result = []
    for item in records:
        values = [item.get(key) for key in keys]
        values_exists = [i not in [None, ''] for i in values]
        exist = values_exists[0]
        if not exist:
            continue
        for e in values_exists[1:]:
            if exist and not e:
                result.append(item)
                break

    return result


def load_abdab(folder):
    file_path = folder / ABDAB_FILE
    result = []
    for item in load_csv(file_path):
        if item['Name'] in SWITCH_MAB:
            item['Name'] = SWITCH_MAB[item['Name']]
        result.append(item)

    return result


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
        heave_v_gene = item['Heavy V Gene']
        if heave_v_gene in ['N/A', 'ND', 'Expected']:
            heave_v_gene = ''
        elif heave_v_gene and '(' in heave_v_gene:
            heave_v_gene = heave_v_gene[:heave_v_gene.index('(')]
        heave_j_gene = item['Heavy J Gene']
        if heave_j_gene in ['N/A', 'ND', 'Expected']:
            heave_j_gene = ''
        elif heave_j_gene and '(' in heave_j_gene:
            heave_j_gene = heave_j_gene[:heave_j_gene.index('(')]
        light_v_gene = item['Light V Gene']
        if light_v_gene in ['N/A', 'ND', 'Expected']:
            light_v_gene = ''
        elif light_v_gene and '(' in light_v_gene:
            light_v_gene = light_v_gene[:light_v_gene.index('(')]
        light_j_gene = item['Light J Gene']
        if light_j_gene in ['N/A', 'ND', 'Expected']:
            light_j_gene = ''
        elif light_j_gene and '(' in light_j_gene:
            light_j_gene = light_j_gene[:light_j_gene.index('(')]

        sources = item['Sources']
        if sources:
            sources = sources[:sources.find('(')]

        if name in mab2Seq:
            print('Duplicated name:', name)
        mab2Seq[name] = {
            'Author': sources,
            'Mab': name,
            'Heavy V Gene': heave_v_gene,
            'VH': vh,
            'VL': vl,
            'V region(H)': '',
            'CDRH3': cdrh3,
            'Heavy J Gene': heave_j_gene,
            'J region(H)': '',
            'Light V Gene': light_v_gene,
            'V region(L)': '',
            'CDRL3': cdrl3,
            'Light J Gene': light_j_gene,
            'J region(L)': '',
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
        mab_list = item['mAb (live EC50, ng/ml or PV: ng/ml)']
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
        mab_list = item['mAb (live EC50, ng/ml or PV: ng/ml)']
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
        if not re.search(r'\d+', refIDs):
            print('RefID dont have publish year', refIDs)

    dup_mab_names = get_dup_key(structures, 'mAb Names')
    print('Duplicated Mab', dup_mab_names)

    dup_pdb = get_dup_key(structures, 'PDB')
    print('Duplicated PDB', dup_pdb)

    structures = [s for s in structures if s['mAb Target'] == 'RBD']

    invalid_records = get_invalid_null_records(structures, [
        'mAb Target', 'Epitope(<=4.5Å)',
        'ACE2-Competing', 'PDB Spike'])

    for item in invalid_records:
        print('Invalid structure record', item)


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


def get_mab_has_problem(structures):
    mab_names = set()
    for item in structures:
        name = item['mAb Names']
        name = name.strip()
        if ';' in name:
            continue

        comment = item['Comment']
        release_date = item['Release date']
        if comment or not release_date:
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
        if not re.search(r'\d+', refIDs):
            print('RefID dont have publish year', refIDs)

    dup_mab = get_dup_key(sequences, 'Mab')
    if dup_mab:
        print('Duplicated mab sequences', dup_mab)

    invalid_hl = get_depend_invalid_null_records(sequences, [
            'VH or VHH', 'IGHV (% Somatic hypermutation)', 'Species',
            'V region(H)', 'CDRH3', 'IGHJ', 'J region(H)'
        ]
    )
    for record in invalid_hl:
        print('Invalid H region', record['Mab'])


def export_sequences_for_alignment(folder, sequences):
    results = []
    for record in sequences:
        new_record = {}
        new_record['Author'] = record['Author']
        new_record['Name'] = record['Mab']
        new_record['Species'] = record['Species']
        new_record['VH or VHH'] = record['VH or VHH']
        new_record['VL'] = record['VL']

        if not new_record['VH or VHH'] and not new_record['VL']:
            continue

        ighv = record.get('IGHV (% Somatic hypermutation)', '')
        ighv = ighv.split('(', 1)[0].strip()
        new_record['IGHV(gene allele)'] = ighv
        ighv = ighv.split('*', 1)[0].strip()
        new_record['Heavy V Gene'] = ighv

        ighj = record.get('IGHJ', '')
        ighj = ighj.split('*', 1)[0].strip()
        new_record['Heavy J Gene'] = ighj

        new_record['IGHJ(gene allele)'] = ''

        iglv = record.get('IGLV (% Somatic hypermutation)', '')
        iglv = iglv.split('(', 1)[0].strip()
        new_record['IGLV(gene allele)'] = iglv
        iglv = iglv.split('*', 1)[0].strip()
        new_record['Light V Gene'] = iglv

        iglj = record.get('IGLJ', '')
        iglj = iglj.split('*', 1)[0].strip()
        new_record['Light J Gene'] = iglj

        new_record['IGLJ(gene allele)'] = ''

        cdrh3 = record.get('CDRH3', '')
        cdrh3 = cdrh3.split('(', 1)[0].strip()
        new_record['CDRH3'] = cdrh3

        cdrl3 = record.get('CDRL3', '')
        cdrl3 = cdrl3.split('(', 1)[0].strip()
        new_record['CDRL3'] = cdrl3

        results.append(new_record)

    dump_csv(folder / 'Mab-all.csv', results)


def show_mab_diff(list1, list1_name, list2, list2_name):
    for i, name in enumerate(list1):
        if name in RENAME_MAB:
            list1[i] = RENAME_MAB[name]

    for i, name in enumerate(list2):
        if name in RENAME_MAB:
            list2[i] = RENAME_MAB[name]

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
    mab_with_problem = get_mab_has_problem(structures)

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
        name = RENAME_MAB_ABDAB2MAb.get(name, name)
        if name not in mab_with_epitope and name not in mab_with_problem:
            miss_epitope.add(name)
    print('PDB without epitope process:', miss_epitope)


def get_missing_sequence(sequences, mab2Seq):
    seq_mabs = [i['Mab'] for i in sequences]
    seq_mabs = [i if i not in RENAME_MAB_ABDAB2MAb
                else RENAME_MAB_ABDAB2MAb[i] for i in seq_mabs]

    missing_seq = []
    for mab, value in mab2Seq.items():
        if mab not in seq_mabs:
            value['Mab'] = mab
            missing_seq.append(value)
            print('Sequence not in table:', mab)

    return missing_seq


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
            print('Seq Error VH:', name)
            print(base_line['VH'])
            # print(vh)
            # print(base_line['VH'] in vh)
        if vl != base_line['VL']:
            print('Seq Error VL:', name)
            print(base_line['VL'])
            # print(vl)
            # print(base_line['VL'] in vl)
        if base_line['CDRH3'] not in cdrh3:
            print('Seq Error CDRH3:', name)
            print(base_line['CDRH3'])
        if base_line['CDRL3'] not in cdrl3:
            print('Seq Error CDRL3:', name)
            print(base_line['CDRL3'])


@click.command()
@click.argument('folder_path', type=click.Path(exists=True))
def work(folder_path):
    folder = Path(folder_path).absolute()

    mabs = load_csv(folder / 'MAbs.csv')
    structures = load_csv(folder / 'Structures.csv')
    sequences = load_csv(folder / 'Sequences.csv')

    print('Verify structures')
    validate_structures(structures)
    print('#' * 80)
    print('Verify sequences')
    validate_sequences(sequences)
    print('#' * 80)
    export_sequences_for_alignment(folder, sequences)

    mab_names1 = get_mab_from_mabs(mabs)
    mab_names2 = get_mab_from_structures(structures)
    mab_names3 = get_mab_from_sequences(sequences)
    show_mab_diff(mab_names1, 'Mabs.csv', mab_names2, 'Structures.csv')
    print('Should be 0.')
    print('#' * 80)
    show_mab_diff(mab_names1, 'Mabs.csv', mab_names3, 'Sequences.csv')
    print('Should be 0.')
    print('#' * 80)

    mab2RefID_base = get_mab2refIDs_from_mabs(mabs)
    mab2RefID_struct = get_mab2refIDs_from_structures(structures)
    mab2RefID_seq = get_mab2refIDs_from_sequences(sequences)

    print('Structure RefIDs')
    show_error_refID(mab2RefID_base, mab2RefID_struct)
    print('Should be 0.')
    print('#' * 80)

    print('Sequence RefIDs')
    show_error_refID(mab2RefID_base, mab2RefID_seq)
    print('Should be 0.')
    print('#' * 80)

    mab_names2 = get_mab_has_epitope(structures)
    show_mab_diff(mab_names2, 'Structures.csv', mab_names3, 'Sequences.csv')
    print('#' * 80)

    abdab = load_abdab(folder)
    mab2PDB = get_mab2PDB_from_abdab(abdab)
    show_missing_structure(structures, mab2PDB)
    print('#' * 80)

    mab2Seq = get_mab2Seq_from_abdab(abdab)
    show_error_sequences(sequences, mab2Seq)
    # missing_seqs = get_missing_sequence(sequences, mab2Seq)
    # save_file = folder / 'Missing-sequence.csv'
    # dump_csv(save_file, missing_seqs)


if __name__ == '__main__':
    work()

# TODO
# check S trimer PDB file, want to replace RBD PDB file
