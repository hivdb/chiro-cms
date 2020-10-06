import os
import json

from .func_mutannots import yield_mutannots_json


def build_mutannots_json(resource_dir, buildres_dir, **kw):
    for resname, payload, _ in yield_mutannots_json(resource_dir):
        dest_json = os.path.join(
            buildres_dir, 'mutannot-{}.json'.format(resname))
        with open(dest_json, 'w') as fp:
            json.dump(payload, fp, indent=2)
            print('create: {}'.format(dest_json))
