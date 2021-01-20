import os
import json
import click
from pprint import pprint
from collections import defaultdict

BASE_DIR = os.path.dirname(
    os.path.dirname(__file__)
)

CACHE_DIR = os.path.join(BASE_DIR, 'local', 'cache-covidcg-aapcnt')

RESULTS_DIR = os.path.join(BASE_DIR, 'resources', 'covidcg-aapcnt')

SUPPORT_GENES = ['S']


def process_gene(content, gene):
    gene_aa_idx_map = {}
    for k, v in content['metadata_map']['gene_aa_snp'].items():
        if k.split('|')[0] == gene:
            gene_aa_idx_map[v] = k

    gene_mut_counts = defaultdict(int)

    total = len(content['case_data'])

    for i in content['case_data']:
        for idx in i['gene_aa_snp_str']:
            mut = gene_aa_idx_map.get(idx)
            if not mut:
                continue
            gene_mut_counts[mut] += 1

    result = []
    for k, v in gene_mut_counts.items():
        gene, pos, ref, mut = k.split('|')
        result.append({
            'gene': gene,
            'position': pos,
            'aa': mut,
            'percent': v / total,
            'count': v,
            'total': total
        })
    result.sort(key=lambda x: int(x['position']))
    return result


@click.command()
@click.argument('genes', nargs=-1, type=click.Choice(SUPPORT_GENES))
def main(genes):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(os.path.join(CACHE_DIR, 'data_package.json')) as fp:
        content = json.load(fp)
        print(content['data_date'])

        for gene in genes:
            result = process_gene(content, gene)
            destfile = os.path.join(RESULTS_DIR, 'aapcnt-{}.json'.format(gene))
            with open(destfile, 'w') as fp:
                json.dump(result, fp, indent=2)


if __name__ == '__main__':
    main()
