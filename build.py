#! /usr/bin/env python

import os
import re
import json
import copy
import subprocess
from collections.abc import Mapping
from importlib import import_module
from datetime import datetime

import ruamel.yaml
from distutils.dir_util import copy_tree, remove_tree

BASEDIR = os.path.dirname(os.path.abspath(__file__))
PAGEDIR = os.path.join(BASEDIR, 'pages')
RESOURCEDIR = os.path.join(BASEDIR, 'resources')
IMAGEDIR = os.path.join(BASEDIR, 'images')
DOWNLOADDIR = os.path.join(BASEDIR, 'downloads')
BUILDDIR = os.path.join(BASEDIR, 'build')
BUILDRESDIR = os.path.join(BUILDDIR, '.resources')
PLUGINDIR = os.path.join(BASEDIR, 'build_plugins')

YAML_PATTERN = re.compile(r'^(.*)\.ya?ml$')

yaml = ruamel.yaml.YAML()


def getmtime(resource_path):
    print(resource_path)
    proc = subprocess.run(
        ['git', 'status', '-s', resource_path],
        stdout=subprocess.PIPE, check=True
    )
    if not proc.stdout:
        proc = subprocess.run(
            ['git', 'log', '-1', '--date', 'unix', resource_path],
            stdout=subprocess.PIPE, check=False
        )
        if proc.returncode == 0:
            for row in proc.stdout.splitlines():
                if row.startswith(b'Date:'):
                    return float(row[5:].strip())
    return os.path.getmtime(resource_path)


def load_resources(data):
    mtime = 0
    data = copy.copy(data)
    if isinstance(data, dict):
        if '_resource' in data:
            rpath = data.pop('_resource')
            resource_path = os.path.join(BUILDRESDIR, rpath)
            if not os.path.isfile(resource_path):
                resource_path = os.path.join(RESOURCEDIR, rpath)
            mtime = max(getmtime(resource_path), mtime)
            with open(resource_path) as fp:
                if resource_path.endswith('.json'):
                    data = json.load(fp)
                elif YAML_PATTERN.match(resource_path):
                    data = yaml.load(fp)
                else:
                    data = fp.read()
        else:
            for key, val in list(data.items()):
                val, nested_mtime = load_resources(val)
                if key == '.' or key.startswith('.<<'):
                    data.pop(key)
                    data.update(val)
                else:
                    data[key] = val
                mtime = max(nested_mtime, mtime)
    elif isinstance(data, list):
        new_data = []
        for val in data:
            val, nested_mtime = load_resources(val)
            new_data.append(val)
            mtime = max(nested_mtime, mtime)
        data = new_data
    return data, mtime


def run_plugins():
    for pyfile in os.listdir(PLUGINDIR):
        if not os.path.isfile(os.path.join(PLUGINDIR, pyfile)) or \
                not re.search(r'^build_.+\.py$', pyfile):
            continue
        modulename, _ = os.path.splitext(pyfile)
        module = import_module('build_plugins.{}'.format(modulename))
        getattr(module, modulename)(
            base_dir=BASEDIR,
            page_dir=PAGEDIR,
            resource_dir=RESOURCEDIR,
            image_dir=IMAGEDIR,
            download_dir=DOWNLOADDIR,
            build_dir=BUILDDIR,
            buildres_dir=BUILDRESDIR,
            plugin_dir=PLUGINDIR
        )


def main():
    os.makedirs(BUILDDIR, exist_ok=True)
    os.makedirs(BUILDRESDIR, exist_ok=True)
    run_plugins()
    for folder, _, files in os.walk(PAGEDIR):
        for filename in files:
            if not YAML_PATTERN.match(filename):
                continue
            yamlpath = os.path.join(folder, filename)
            mtime = getmtime(yamlpath)
            rel_yamlpath = os.path.relpath(yamlpath, BASEDIR)
            jsonpath = os.path.join(
                BUILDDIR, YAML_PATTERN.sub(r'\1.json', rel_yamlpath))
            os.makedirs(os.path.dirname(jsonpath), exist_ok=True)
            with open(yamlpath) as yamlfp, open(jsonpath, 'w') as jsonfp:
                data = yaml.load(yamlfp)
                data, new_mtime = load_resources(data)
                mtime = max(new_mtime, mtime)
                if isinstance(data, Mapping):
                    data['lastModified'] = (
                        datetime.utcfromtimestamp(mtime).isoformat() + 'Z'
                    )
                json.dump(data, jsonfp)
                print('create: {}'.format(jsonpath))
    copy_tree(IMAGEDIR, os.path.join(BUILDDIR, 'images'))
    copy_tree(DOWNLOADDIR, os.path.join(BUILDDIR, 'downloads'))
    remove_tree(BUILDRESDIR)


if __name__ == '__main__':
    main()
