#! /usr/bin/env python

import csv
import json
import click
from datetime import datetime
# import textwrap
# import ruamel.yaml
# from ruamel.yaml.comments import CommentedSeq
# from ruamel.yaml.scalarstring import LiteralScalarString


@click.command()
@click.argument('input_file', type=click.File())
@click.argument('output_file', type=click.File('w'))
@click.argument('version', type=str)
def mutcmt2json(input_file, output_file, version):
    rows = []  # CommentedSeq()
    for idx, row in enumerate(csv.DictReader(input_file)):
        row['position'] = int(row['position'])
        last_update = datetime.strptime(
            row['date_updated'] or row['date_added'],
            '%Y-%m-%d'
        ).strftime('%b %-d, %Y')
        row['comment'] = (
            '{} ({})'
            .format(row['comment'], last_update)
            .replace('  ', ' ')
        )
        # row['comments'] = LiteralScalarString(
        #     '\n'.join(textwrap.wrap(
        #         row['comments'],
        #         width=75,
        #         break_on_hyphens=False,
        #         break_long_words=False
        #     ))
        # )
        rows.append(row)
        # if idx > 0:
        #     rows.yaml_set_comment_before_after_key(idx, "\n\n", 2)
    # yaml = ruamel.yaml.YAML()
    # yaml.dump(rows, output_file)
    json.dump({
        'version': version,
        'payload': rows
    }, output_file)


if __name__ == '__main__':
    mutcmt2json()
