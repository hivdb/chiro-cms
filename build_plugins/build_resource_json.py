import os
import re
import json
import ruamel.yaml


def build_resource_json(resource_dir, buildres_dir, **kw):
    for resyaml in os.listdir(resource_dir):
        resyaml_path = os.path.join(resource_dir, resyaml)
        if not os.path.isfile(resyaml_path) or \
                not re.search(r'.+(?<!-data)\.ya?ml$', resyaml):
            continue
        with open(resyaml_path, encoding='utf-8-sig') as fp:
            res_config = ruamel.yaml.load(fp, Loader=ruamel.yaml.Loader)
        resname, _ = os.path.splitext(resyaml)
        data_dir = '{}-data'.format(resname)
        data_dir = os.path.join(resource_dir, data_dir)
        data = []
        data_loaded = False
        for suffix in ('yml', 'yaml', 'json'):
            yaml_path = os.path.extsep.join([data_dir, suffix])
            if os.path.isfile(yaml_path):
                data_loaded = True
                with open(yaml_path, encoding='utf-8-sig') as fp:
                    data = ruamel.yaml.load(fp, Loader=ruamel.yaml.Loader)
                    break
        if not data_loaded:
            for basedir, _, yaml_paths in os.walk(data_dir):
                for yaml_path in yaml_paths:
                    yaml_path = os.path.join(basedir, yaml_path)
                    if not os.path.isfile(yaml_path) or \
                            not re.search(r'\.ya?ml$', yaml_path):
                        continue
                    with open(yaml_path, encoding='utf-8-sig') as fp:
                        data_loaded = True
                        group = ruamel.yaml.load(fp, Loader=ruamel.yaml.Loader)
                        data.append(group)
        dest_json = os.path.join(buildres_dir, '{}.json'.format(resname))
        with open(dest_json, 'w') as fp:
            if data_loaded:
                res_config['data'] = data
            json.dump(res_config, fp, indent=2)
            print('create: {}'.format(dest_json))
