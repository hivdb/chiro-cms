import os
import json
from itertools import groupby

from .func_mutannots import (get_positions,
                             get_amino_acids,
                             load_config_and_data)


def build_mutannots_table(resource_dir, buildres_dir, **kw):
    for resname, geneconfig, annotdata_lookup in \
            load_config_and_data(resource_dir):

        refseq = geneconfig['refSequence']

        coldefs = [
            {'name': 'position',
             'label': 'Pos',
             'sort': 'numeric'},
            {'name': 'refAA',
             'label': 'Ref'}
        ]
        rows = {
            pos0 + 1: {
                'position': pos0 + 1,
                'refAA': aa
            } for pos0, aa in enumerate(refseq)
        }
        for cat in geneconfig['annotCategories']:
            cat_annots = cat.pop('annotNames')
            for annot_name in cat_annots:
                annotdata = annotdata_lookup[annot_name]
                annot_label = annotdata.get(
                    'tableColLabel',
                    annotdata.get('label', annot_name)
                )
                level = annotdata['level']

                pos_by_cite = [
                    set(get_positions(cite)) for cite in annotdata['citations']
                ]

                if all(pos == pos_by_cite[0] for pos in pos_by_cite):
                    pos_by_cite = pos_by_cite[:1]
                    col_names = [annot_label]
                    coldefs.append({
                        'name': annot_label,
                        'label': annot_label,
                        'none': ''
                    })
                else:
                    col_names = [
                        '{}-{author}-{year}'.format(annot_label, **cite)
                        for cite in annotdata['citations']
                    ]
                    coldefs.extend([{
                        'name': ('{}-{author}-{year}'
                                 .format(annot_label, **cite)),
                        'label': ('{}\n({}{})'
                                  .format(annot_label,
                                          cite['author'],
                                          cite['year'] - 2000)),
                        'none': ''
                    } for cite in annotdata['citations']])

                for col_name, posset in zip(col_names, pos_by_cite):
                    if level == 'position':
                        for subgroup in annotdata['positions']:
                            for pos in get_positions(subgroup):
                                if pos not in posset:
                                    continue
                                row = rows[pos]
                                row[col_name] = subgroup['subgroup']
                    else:
                        for citedata in annotdata['citations']:
                            for pos, posaas in groupby(
                                get_amino_acids(citedata),
                                lambda posaa: posaa[0]
                            ):
                                if pos not in posset:
                                    continue
                                row = rows[pos]
                                aas = [aa for _, aa in posaas]
                                row[col_name] = '/'.join(sorted(aas))

        rows = list(rows.values())
        dest_json = os.path.join(
            buildres_dir, 'mutannot-{}-table.json'.format(resname))
        with open(dest_json, 'w') as fp:
            json.dump({
                'columnDefs': coldefs,
                'data': rows
            }, fp, indent=2)
            print('create: {}'.format(dest_json))
