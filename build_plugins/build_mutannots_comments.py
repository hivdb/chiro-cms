import re
import os
import json

from .func_mutannots import yield_mutannots_json


def parse_mutations(mutations_text):
    for pos, aas in re.findall(
        r'(\d+)([A-Z]+|ins|del|stop)(?:$|,)',
        mutations_text
    ):
        aas = (aas
               .replace('ins', 'i')
               .replace('del', 'd')
               .replace('stop', '*'))
        yield int(pos), aas


def build_mutannots_comments(resource_dir, download_dir, buildres_dir, **kw):
    with open(os.path.join(
        download_dir,
        'mutation-comments',
        'latest.json'
    )) as fp:
        all_comments = json.load(fp)

    for resname, payload, geneconfig in yield_mutannots_json(resource_dir):
        gene = payload['gene']
        gene_comments = {}
        for comment in all_comments['payload']:
            if comment['gene'] != gene:
                continue
            for pos, aas in parse_mutations(comment['mutations']):
                gene_comments.setdefault(pos, []).append(comment['comment'])
        gene_comments = [{
            'position': pos,
            'comment': '\n'.join(cmts)
        } for pos, cmts in gene_comments.items()]
        dest_json = os.path.join(
            buildres_dir, '{}-comments.json'.format(resname))
        with open(dest_json, 'w') as fp:
            json.dump({
                'columnDefs': [{
                    'name': 'position',
                    'label': 'Pos',
                    'sort': 'numeric'
                },  {
                    # 'name': 'triggeredAAs',
                    # 'label': 'AAs'
                    # }, {
                    'name': 'comment',
                    'label': 'Comment',
                    'textAlign': 'left',
                    'sortable': False
                }],
                'data': gene_comments,
                'references': ''
            }, fp, indent=2)
            print('create: {}'.format(dest_json))
