import click
import csv
import yaml
from pathlib import Path


def trim_records(records):
    for record in records:
        for k, v in record.items():
            record[k] = v.strip()

    return records


def load_csv(file_path, trim=True):
    records = []
    with open(file_path, encoding='utf-8-sig') as fd:
        for record in csv.DictReader(fd):
            records.append(record)

    if trim:
        return trim_records(records)
    else:
        return records


def load_yaml(file_path):
    with file_path.open() as fd:
        return yaml.load(fd, Loader=yaml.Loader)


def dump_yaml(file_path, yaml_obj):
    with file_path.open('w', encoding='utf-8') as fd:
        yaml.dump(yaml_obj, fd, allow_unicode=True)


def check_mab_no_seq(se, st):
    mabs = []
    for item in se:
        mab = item['Mab']
        mabs.append(mab)

    no_seq = []
    for item in st:
        mab = item['mAb Names']
        if not mab:
            continue
        if mab not in mabs:
            no_seq.append(mab)

    for i in no_seq:
        print(i, 'noseq')


@click.command()
@click.argument('folder_path', type=click.Path(exists=True))
@click.argument('output_path',
                type=click.Path(dir_okay=True, file_okay=False))
def import_antibody_metainfo(folder_path, output_path):
    folder = Path(folder_path).absolute()

    st = load_csv(folder / 'Structures.csv')
    se = load_csv(folder / 'Sequences.csv')

    se_map = {}
    for item in se:
        mab = item['Mab']
        IGHV = item['IGHV (% Somatic hypermutation)']
        if not IGHV:
            continue
        IGHV = IGHV.split()[0].split('*')[0].replace('IGHV', '')
        se_map[mab] = IGHV

    result = []
    for item in st:
        mab = item['mAb Names']
        PDB = item['PDB']
        if not PDB:
            continue
        target = item['mAb Target']
        if target not in ['RBD', 'NTD']:
            continue

        ace2_competing = item['ACE2-Competing']
        ace2_overlap = item['%ACE2 Overlap']
        s_trim = item['PDB Spike']
        if not ace2_competing:
            continue
        if not str(ace2_overlap):
            continue
        if not s_trim:
            continue
        result.append({
            'mAb': mab,
            'PDB': PDB,
            'PDB_s_trimer': ('√' if s_trim == 'S Trimer' else ''),
            'target': target,
            'IGHV': se_map[mab],
            'ace2_competing': ('√' if ace2_competing == 'Yes' else ''),
            'ace2_overlap': int(ace2_overlap)
        })

    save_path = Path(output_path) / 'antibody-meta.yml'
    dump_yaml(save_path, result)

    # check_mab_no_seq(se, st)


if __name__ == '__main__':
    import_antibody_metainfo()
