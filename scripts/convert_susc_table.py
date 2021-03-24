#! /usr/bin/env python
import click
import json


@click.command()
@click.argument('input_file', type=click.Path())
@click.argument('output_file', type=click.Path())
def convert_susc_table(input_file, output_file):
    with open(input_file, encoding='utf-8') as fd:
        data = json.load(fd)
    result = []
    for item in data:
        if not item.get('assays'):
            result.append(item)
            continue

        new_assays = []
        for rec in item.get('assays'):
            refname = rec['reference']
            refID = refname.replace('*', '').replace('†', '')
            rec['reference'] = "{}[^{}]".format(refname, refID)
            new_assays.append(rec)
        item['assays'] = new_assays
        result.append(item)

    with open(output_file, 'w', encoding='utf-8') as fd:
        json.dump(result, fd, indent=4, ensure_ascii=False)


if __name__ == '__main__':
    convert_susc_table()
