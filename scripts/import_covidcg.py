import os
import json
import click  # type: ignore
import ruamel.yaml  # type: ignore
from collections import defaultdict

yaml = ruamel.yaml.YAML()

BASE_DIR = os.path.dirname(
    os.path.dirname(__file__)
)

CACHE_DIR = os.path.join(BASE_DIR, 'local', 'cache-covidcg-aapcnt')

RESULTS_DIR = os.path.join(BASE_DIR, 'resources', 'covidcg-aapcnt')

SUPPORT_GENES = ['S']

GENEREF_S = (
    'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGT'
    'NGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYH'
    'KNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQG'
    'FSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCA'
    'LDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADY'
    'SVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNS'
    'NNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVV'
    'VLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEI'
    'LDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVN'
    'NSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVS'
    'MTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFS'
    'QILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSA'
    'LLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQ'
    'DVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASA'
    'NLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVF'
    'VSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGD'
    'ISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCC'
    'SCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*'
)

GENEREFS = {
    'S': GENEREF_S
}

GENE_MUTANNOTS = {
    'S': [{
        'cutoff': 0.001,
        'path': os.path.join(
            BASE_DIR, 'resources',
            'mutannots', 'spike', 'var1000.yml')
    }, {
        'cutoff': 0.0005,
        'path': os.path.join(
            BASE_DIR, 'resources',
            'mutannots', 'spike', 'var2000.yml')
    }]
}


def inline_list(*lst):
    ret = ruamel.yaml.comments.CommentedSeq(lst)
    ret.fa.set_flow_style()
    return ret


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
        mut = mut.replace('-', 'd')
        result.append({
            'gene': gene,
            'position': int(pos),
            'aa': mut,
            'percent': v / total,
            'count': v,
            'total': total
        })
    result.sort(key=lambda x: int(x['position']))
    return result


def update_mutannots(gene, aapcnt):
    generef = GENEREFS[gene]
    for config in GENE_MUTANNOTS[gene]:
        cutoff = config['cutoff']
        path = config['path']
        if os.path.isfile(path):
            with open(path) as fp:
                data = yaml.load(fp)
                if not data:
                    raise RuntimeError('File {} is empty'.format(path))
        else:
            data = {
                'name': 'var{}'.format(round(1 / cutoff)),
                'label': 'COVID CG≥{:g}%'.format(cutoff * 100),
                'level': 'aminoAcid',
                'hideCitations': True,
                'colorRules': [],
                'citations': [{
                    'doi': '10.2807/1560-7917.ES.2017.22.13.30494',
                    'author': 'Shu',
                    'year': 2017,
                    'section': 'https://covidcg.org',
                    'aminoAcids': []
                }]
            }
        posaas = defaultdict(list)
        pcnt_pairs = {}
        for one in aapcnt:
            pos = one['position']
            aa = one['aa']
            if one['percent'] < cutoff:
                continue
            if aa == generef[pos - 1]:
                # skip reference AA
                continue
            posaas[pos].append(aa)
            pcnt_pairs['{}{}'.format(pos, aa)] = one['percent']
        posaas = ['{}{}'.format(pos, ''.join(aas))
                  for pos, aas in posaas.items()]
        data['citations'][0]['aminoAcids'] = inline_list(*posaas)
        data['citations'][0]['aminoAcidAttrs'] = [{
            'attr': 'prevalence',
            'pairs': pcnt_pairs
        }]
        with open(path, 'w') as fp:
            yaml.dump(data, fp)
            print('Updated {}'.format(path))


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
            update_mutannots(gene, result)


if __name__ == '__main__':
    main()
