import os
import re
import json
import click

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
    return [round(freq, 4) for freq in freqs]


def iter_percents(gene, branch_freqs, branch_lookup, pivots):
    num_pivots = len(pivots)
    all_branches = set(branch_freqs.keys())
    for pos0, refaa in enumerate(GENEREFS[gene]):
        pos = pos0 + 1
        lookup = branch_lookup.get((gene, pos), {})
        ref_branches = set(all_branches)
        for aa, var_branches in lookup.items():
            if aa == refaa:
                continue
            ref_branches -= var_branches
            freqs = get_freqs(var_branches, branch_freqs, num_pivots)
            yield {
                'gene': gene,
                'position': pos,
                'aa': aa,
                'percents': list(zip(pivots, freqs))
            }
        freqs = get_freqs(ref_branches, branch_freqs, num_pivots)
        yield {
          'gene': gene,
          'position': pos,
          'aa': refaa,
          'percents': list(zip(pivots, freqs))
        }


@click.command()
@click.argument('genes', nargs=-1, type=click.Choice(SUPPORT_GENES))
def main(genes):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(os.path.join(CACHE_DIR, 'global-frequencies.json')) as fp:
        branch_freqs = json.load(fp)
        pivots = branch_freqs.pop('pivots')
        branch_freqs.pop('generated_by')
    branches = set(branch_freqs.keys())
    with open(os.path.join(CACHE_DIR, 'global.json')) as fp:
        data = json.load(fp)
        branch_lookup = build_branch_lookup(data['tree'], branches)
    for gene in genes:
        destfile = os.path.join(RESULTS_DIR, 'aapcnt-{}.json'.format(gene))
        with open(destfile, 'w') as fp:
            payload = list(iter_percents(gene, branch_freqs,
                                         branch_lookup, pivots))
            json.dump(payload, fp, indent=2)
            click.echo('Created {}'.format(destfile))


if __name__ == '__main__':
    main()
