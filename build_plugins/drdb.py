import yaml
import sqlite3
import requests
from tempfile import NamedTemporaryFile
from contextlib import contextmanager

from typing import Generator

from pathlib import Path


BASE_DIR = Path(__file__).parents[1]
DRDB_VERSION_YAML = BASE_DIR / 'resources' / 'drdb-version.yml'


@contextmanager
def open_drdb(
    target: str = 'full'
) -> Generator[sqlite3.Connection, None, None]:
    with DRDB_VERSION_YAML.open(encoding='UTF-8-sig') as fp:
        drdb_version_yaml = yaml.load(fp, Loader=yaml.Loader)
    drdb_version: str = drdb_version_yaml[target]
    resp = requests.get(
        'https://s3-us-west-2.amazonaws.com/cms.hivdb.org'
        '/covid-drdb/covid-drdb-{}.db'.format(drdb_version)
    )
    with NamedTemporaryFile() as dbfile:
        dbfile.write(resp.content)
        dbfile.seek(0)
        yield sqlite3.connect(dbfile.name)
