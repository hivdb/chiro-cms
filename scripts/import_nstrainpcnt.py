import os
import re
import json
import click
from collections import Counter
from copy import deepcopy
from pprint import pprint
from datetime import date

BASE_DIR = os.path.dirname(
    os.path.dirname(__file__)
)

CACHE_DIR = os.path.join(BASE_DIR, 'local', 'cache-nstrain-aapcnt')

RESULTS_DIR = os.path.join(BASE_DIR, 'resources', 'nstrain-aapcnt')

SUPPORT_GENES = ['RdRP', 'S']

GENEREF_RDRP = (
    'SADAQSFLNRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKTNCCRFQEKDEDDNLIDSYFVVK'
    'RHTFSNYQHEETIYNLLKDCPAVAKHDFFKFRIDGDMVPHISRQRLTKYTMADLVYALRHFDEGNCDTLKEIL'
    'VTYNCCDDDYFNKKDWYDFVENPDILRVYANLGERVRQALLKTVQFCDAMRNAGIVGVLTLDNQDLNGNWYDF'
    'GDFIQTTPGSGVPVVDSYYSLLMPILTLTRALTAESHVDTDLTKPYIKWDLLKYDFTEERLKLFDRYFKYWDQ'
    'TYHPNCVNCLDDRCILHCANFNVLFSTVFPPTSFGPLVRKIFVDGVPFVVSTGYHFRELGVVHNQDVNLHSSR'
    'LSFKELLVYAADPAMHAASGNLLLDKRTTCFSVAALTNNVAFQTVKPGNFNKDFYDFAVSKGFFKEGSSVELK'
    'HFFFAQDGNAAISDYDYYRYNLPTMCDIRQLLFVVEVVDKYFDCYDGGCINANQVIVNNLDKSAGFPFNKWGK'
    'ARLYYDSMSYEDQDALFAYTKRNVIPTITQMNLKYAISAKNRARTVAGVSICSTMTNRQFHQKLLKSIAATRG'
    'ATVVIGTSKFYGGWHNMLKTVYSDVENPHLMGWDYPKCDRAMPNMLRIMASLVLARKHTTCCSLSHRFYRLAN'
    'ECAQVLSEMVMCGGSLYVKPGGTSSGDATTAYANSVFNICQAVTANVNALLSTDGNKIADKYVRNLQHRLYEC'
    'LYRNRDVDTDFVNEFYAYLRKHFSMMILSDDAVVCFNSTYASQGLVASIKNFKSVLYYQNNVFMSEAKCWTET'
    'DLTKGPHEFCSQHTMLVKQGDDYVYLPYPDPSRILGAGCFVDDIVKTDGTLMIERFVSLAIDAYPLTKHPNQE'
    'YADVFHLYLQYIRKLHDELTGHMLDMYSVMLTNDNTSRYWEPEFYEAMYTPHTVLQ'
)

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
    'RdRP': GENEREF_RDRP,
    'S': GENEREF_S
}


def parse_mutation(gene, mut):
    match = re.match(r'^[A-Z*-](\d+)([A-Z*-])$', mut)
    if not match:
        click.echo('Unable to recongnize mutation {}'.format(mut))
        return
    pos, aa = match.groups()
    pos = int(pos)
    if gene == 'ORF1a' and pos > 4392:
        return 'RdRP', pos - 4392, aa
    elif gene == 'ORF1b' and pos < 924:
        return 'RdRP', pos + 9, aa
    return gene, pos, aa


def parse_mutations(mutations):
    mutobjs = set()
    for gene, muts in mutations.items():
        if gene == 'nuc':
            continue
        for mut in muts:
            mutobj = parse_mutation(gene, mut)
            if mutobj:
                mutobjs.add(mutobj)
    return mutobjs


def _build_branch_lookup(node, keep_branches, parent_mutations=None):
    branch_attrs = node.get('branch_attrs')
    if branch_attrs and 'mutations' in branch_attrs:
        mutations = parse_mutations(branch_attrs['mutations'])
    else:
        mutations = set()

    if parent_mutations:
        mutations |= parent_mutations
    if 'name' in node and node['name'] in keep_branches:
        yield node['name'], mutations
    if 'children' in node:
        for child_node in node['children']:
            yield from _build_branch_lookup(
                child_node, keep_branches, mutations)


def build_branch_lookup(tree, keep_branches):
    lookup = {}
    for name, mutations in _build_branch_lookup(tree, keep_branches):
        for gene, pos, aa in mutations:
            if (gene, pos) not in lookup:
                lookup[(gene, pos)] = {}
            if aa not in lookup[(gene, pos)]:
                lookup[(gene, pos)][aa] = set()
            lookup[(gene, pos)][aa].add(name)
    return lookup


def get_freqs(branch_names, branch_freqs, num_pivots):
    freqs = [0] * num_pivots
    for branch_name in branch_names:
        for idx, freq in enumerate(branch_freqs[branch_name]['frequencies']):
            freqs[idx] += freq
    return freqs


def iter_percents(gene, branch_freqs, branch_lookup, pivots):
    num_pivots = len(pivots)
    all_branches = set(branch_freqs.keys())
    min_var = 1e-10
    for pos0, refaa in enumerate(GENEREFS[gene]):
        pos = pos0 + 1
        lookup = branch_lookup.get((gene, pos), {})
        ref_branches = set(all_branches)
        total_freqs = [0] * len(pivots)
        posaas = []
        for aa, var_branches in lookup.items():
            if aa == refaa:
                continue
            ref_branches -= var_branches
            freqs = get_freqs(var_branches, branch_freqs, num_pivots)
            total_freqs = [a + b for a, b in zip(total_freqs, freqs)]
            posaas.append({
                'gene': gene,
                'position': pos,
                'aa': aa,
                'percents': freqs
            })
        freqs = get_freqs(ref_branches, branch_freqs, num_pivots)
        total_freqs = [a + b for a, b in zip(total_freqs, freqs)]
        posaas.append({
          'gene': gene,
          'position': pos,
          'aa': refaa,
          'percents': freqs
        })
        len_posaas = len(posaas)
        for one in posaas:
            one['percents'] = list(zip(pivots, [
                round((a + min_var) / (min_var * len_posaas + b), 4)
                for a, b in zip(one['percents'], total_freqs)
            ]))
            yield one


def numeric2date(pivots):
    results = []
    for one in pivots:
        year = int(one)
        offset = one - year
        days = date(year + 1, 1, 1) - date(year, 1, 1)
        results.append((date(year - 1, 12, 31) + days * offset).isoformat())
    return results


# Copy functions from auapice/src/actions/frequencies.js
def generate_tree(subtree):
    result = []
    result.append(subtree)
    for item in subtree.get('children', []):
        result.extend(generate_tree(item))
    return result


def parseMutPos(mut):
    return int(mut[1:-1])


def setGenotype(nodes, prot, positions, refSequence):
    nPositions = len(positions)
    ancState = [None for i in positions]
    ancNodes = [[] for i in positions]

    def _setGenotype_recurse(node, state):
        nonlocal nPositions
        nonlocal ancState
        nonlocal ancNodes
        newState = state
        mutations_data = None
        if (node.get('branch_attrs')
                and node['branch_attrs'].get('mutations')
                and node['branch_attrs']['mutations'].get(prot)):
            mutations_data = node['branch_attrs']['mutations'].get(prot)

        if mutations_data:
            for mut in mutations_data:
                mPos = parseMutPos(mut)
                for j in range(nPositions):
                    if (positions[j] == mPos):
                        if (not ancState[j]):
                            ancState[j] = mut[0]
                        newState[j] = mut[-1]

        node['currentGt'] = deepcopy(newState)
        for j in range(nPositions):
            if not newState[j]:
                ancNodes[j].append(node)

        if node.get('children'):
            for child in node['children']:
                _setGenotype_recurse(child, deepcopy(newState))

    _setGenotype_recurse(nodes[0], [None for i in positions])

    if refSequence:
        for i, inferredValue in enumerate(ancState):
            pos = positions[i]
            _, orig_pos = convert_gene_pos_back(prot, pos)
            ref_AA = refSequence[orig_pos - 1]
            if not inferredValue:
                ancState[i] = ref_AA
            elif inferredValue != ref_AA:
                print(
                    'Mismatch between inferred ancestral state for'
                    f'{prot}@{positions[i]} of ${inferredValue}'
                    f' and ${ref_AA}'
                    )

    for i in range(nPositions):
        for node in ancNodes[i]:
            node['currentGt'][i] = ancState[i]

    for node in nodes:
        node['currentGt'] = ' / '.join(node.get('currentGt', []))

    return nodes


def computeMatrixFromRawData(data, pivots, nodes, categories):
    unassigned_label = "unassigned"
    categories.append(unassigned_label)
    matrix = {}
    pivotsLen = len(pivots)
    for i in categories:
        matrix[i] = [0] * pivotsLen

    debugPivotTotals = [0] * pivotsLen
    for d in data:
        cat = nodes[d['idx']].get('currentGt')
        if not cat:
            cat = unassigned_label
        for i in range(pivotsLen):
            matrix[cat][i] += d['values'][i]
            debugPivotTotals[i] += d['values'][i]

    nCategories = len(matrix.keys())
    minVal = 1e-10
    for cat in matrix.keys():
        for idx, norm in enumerate(debugPivotTotals):
            matrix[cat][idx] = (matrix[cat][idx] + minVal) / (
                nCategories * minVal + norm)

    if sum(matrix[unassigned_label]) == 0:
        del matrix[unassigned_label]

    return matrix


def processFrequenciesJSON(rawJSON, tree):
    subtree = []
    for idx, item in enumerate(tree):
        if item.get('children'):
            continue
        else:
            subtree.append((idx, item))

    result = []
    for idx, item in subtree:
        if not rawJSON.get(item['name']):
            continue
        result.append({
            'idx': idx,
            'values': rawJSON[item['name']]['frequencies'],
            'weight': rawJSON[item['name']].get('weight')
        })

    return result


def convert_gene_pos(gene, pos):
    if gene == 'RdRP':
        if pos < 10:
            return 'ORF1a', pos + 4392
        else:
            return 'ORF1b', pos - 9
    return gene, pos


def convert_gene_pos_back(gene, pos):
    if gene == 'ORF1a' and pos > 4392:
        return 'RdRP', pos - 4392
    elif gene == 'ORF1b' and pos < 924:
        return 'RdRP', pos + 9
    return gene, pos


def process_gene(gene, generef, data, pivots, tree):
    posaas = []
    for offset, refaa in enumerate(generef):
        pos = offset + 1
        cal_gene, calc_pos = convert_gene_pos(gene, pos)
        tree = setGenotype(tree, cal_gene, [calc_pos], generef)

        all_muts = []
        for item in tree:
            all_muts.append(item['currentGt'])

        all_muts = list(set(all_muts))

        matrix = computeMatrixFromRawData(
            data,
            pivots,
            tree,
            all_muts
        )

        for key, value in matrix.items():
            if key == 'unassigned':
                continue
            posaas.append({
                'gene': gene,
                'position': pos,
                'aa': key,
                'percents': list(zip(pivots, value))
                })

    return posaas


@click.command()
@click.argument('genes', nargs=-1, type=click.Choice(SUPPORT_GENES))
def main(genes):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(os.path.join(CACHE_DIR, 'global.json')) as fp:
        data = json.load(fp)
        tree = generate_tree(data['tree'])

    with open(os.path.join(CACHE_DIR, 'global-frequencies.json')) as fp:
        branch_freqs = json.load(fp)
        pivots = numeric2date(branch_freqs.pop('pivots'))

    data = processFrequenciesJSON(branch_freqs, tree)

    for gene in genes:
        generef = GENEREFS[gene]
        posaas = process_gene(gene, generef, data, pivots, tree)

        destfile = os.path.join(RESULTS_DIR, 'aapcnt-{}.json'.format(gene))
        with open(destfile, 'w') as fp:
            json.dump(posaas, fp, indent=2)


    # with open(os.path.join(CACHE_DIR, 'global-frequencies.json')) as fp:
    #     branch_freqs = json.load(fp)
    #     pivots = numeric2date(branch_freqs.pop('pivots'))
    #     branch_freqs.pop('generated_by')
    # branches = set(branch_freqs.keys())
    # with open(os.path.join(CACHE_DIR, 'global.json')) as fp:
    #     data = json.load(fp)
    #     branch_lookup = build_branch_lookup(data['tree'], branches)

    # for gene in genes:
    #     destfile = os.path.join(RESULTS_DIR, 'aapcnt-{}.json'.format(gene))
    #     with open(destfile, 'w') as fp:
    #         payload = list(iter_percents(gene, branch_freqs,
    #                                      branch_lookup, pivots))
    #         json.dump(payload, fp, indent=2)
    #         click.echo('Created {}'.format(destfile))


if __name__ == '__main__':
    main()
