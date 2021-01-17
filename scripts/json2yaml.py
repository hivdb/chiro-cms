#! /usr/bin/env python

import json
import click
import ruamel.yaml


@click.command()
@click.argument('input_file', type=click.File())
@click.argument('output_file', type=click.File('w'))
def json2yaml(input_file, output_file):
    yaml = ruamel.yaml.YAML()
    data = json.load(input_file)
    yaml.dump(data, output_file)


if __name__ == '__main__':
    json2yaml()
