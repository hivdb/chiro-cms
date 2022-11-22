import json
import yaml
import sqlite3
import requests
from tempfile import NamedTemporaryFile
from contextlib import contextmanager
from itertools import groupby

from typing import Generator

from pathlib import Path


BASE_DIR = Path(__file__).parents[1]
SARS2_DRMS_YAML = BASE_DIR / 'pages' / 'sars2-drms.yml'
RESIST_MUTS_JSON = (
    BASE_DIR / 'downloads' /
    'resistance-mutations' / 'latest.json'
)

DRUG_CLASSES = ['MAB', '_3CLPI', 'RdRPI']


@contextmanager
def open_drdb() -> Generator[sqlite3.Connection, None, None]:
    with SARS2_DRMS_YAML.open(encoding='UTF-8-sig') as fp:
        drms_yaml = yaml.load(fp, Loader=yaml.Loader)
    drdb_version: str = drms_yaml['drdbVersion']
    resp = requests.get(
        'https://s3-us-west-2.amazonaws.com/cms.hivdb.org'
        '/covid-drdb/covid-drdb-{}.db'.format(drdb_version)
    )
    with NamedTemporaryFile() as dbfile:
        dbfile.write(resp.content)
        dbfile.seek(0)
        yield sqlite3.connect(dbfile.name)


def main() -> None:
    header = ['gene', 'position', 'aa', 'invivo_selected']

    with open_drdb() as conn:
        cur = conn.cursor()
        cur.execute("""
            SELECT
              gene,
              position,
              amino_acid,
              EXISTS (
                SELECT 1 FROM resistance_mutation_attributes RMA
                WHERE
                  RMA.gene = RM.gene AND
                  RMA.position = RM.position AND
                  RMA.amino_acid = RM.amino_acid AND
                  RMA.col_name = 'INVIVO'
              ) AS invivo_selected
            FROM resistance_mutations RM
            ORDER BY gene, position, amino_acid
        """)
        results = [
            dict(zip(header, row))
            for row in cur.fetchall()
        ]
    payload = {}
    for gene, irows in groupby(results, key=lambda row: row['gene']):
        rows = []
        for row in irows:
            row['invivo_selected'] = bool(row['invivo_selected'])
            rows.append(row)
        if gene == 'S':
            payload['MAB'] = rows
        elif gene == '_3CLpro':
            payload['_3CLPI'] = rows
        elif gene == 'RdRP':
            payload['RdRPI'] = rows
        else:
            raise RuntimeError(
                'Unsupport gene {!r} supplied, this script may need an update'
                .format(gene)
            )
    payload = dict(sorted(
        payload.items(),
        key=lambda item: DRUG_CLASSES.index(item[0])
    ))
    with RESIST_MUTS_JSON.open('w') as fp:
        json.dump(payload, fp, indent=2)
        print(fp.name)


if __name__ == '__main__':
    main()
