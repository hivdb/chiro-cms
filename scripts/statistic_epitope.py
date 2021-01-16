import ruamel.yaml
from pathlib import Path
from collections import Counter
import csv

yaml = ruamel.yaml.YAML()

EPITOPE_PATH = (
    Path(__file__).absolute().parent.parent / 'resources' / 'mutannots'
    / 'spike'
)


def iter_epitope_files(folder):
    for i in folder.glob('epitope*'):
        yield i


def load_epitope_file(fd):
    content = yaml.load(fd)
    positions = content['positions']
    all_positions = [i['positions'] for i in positions]
    all_positions = [j for i in all_positions for j in i]
    all_positions = list(set(all_positions))
    return all_positions


def dump_csv(file_path, records, headers=[]):
    if not records:
        return
    if not headers and records:
        headers = records[0].keys()

    with open(file_path, 'w', encoding='utf-8-sig') as fd:
        writer = csv.DictWriter(fd, fieldnames=headers)
        writer.writeheader()
        writer.writerows(records)


def work():
    all_epitope_pos = []
    for fd in iter_epitope_files(EPITOPE_PATH):
        all_epitope_pos.extend(load_epitope_file(fd))

    result = []
    for pos, count in Counter(all_epitope_pos).items():
        result.append({
            'Position': int(pos),
            'epitope count': count
        })

    result.sort(key=lambda x: x['Position'])
    save_path = EPITOPE_PATH / 'epitope.pos.count.csv'
    dump_csv(save_path, result)


if __name__ == '__main__':
    work()
