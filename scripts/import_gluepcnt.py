import os
import csv
import json
import math
import click
import requests
from collections import Counter

from diskcache import FanoutCache

BASE_DIR = os.path.dirname(
    os.path.dirname(__file__)
)

CACHE_DIR = os.path.join(BASE_DIR, 'local', 'cache-glue-aapcnt')

RESULTS_DIR = os.path.join(BASE_DIR, 'resources', 'glue-aapcnt')


CACHE = FanoutCache(CACHE_DIR)

GENEREF_NSP12 = (
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
    'nsp12': GENEREF_NSP12,
    'S': GENEREF_S
}
GENE_OFFSETS = {
    'nsp12': [{
        'between': (13442, 13467),
        'offset': 13441
    }, {
        'between': (13468, 16236),
        'offset': 13440
    }],
    'S': [{
        'between': (21563, 25384),
        'offset': 21562
    }]
}


URL_PREFIX = 'http://cov-glue.cvr.gla.ac.uk/gluetools-ws/'

URL_NUM_SEQUENCES = (
    URL_PREFIX +
    "project/cov/custom-table-row/"
    "cov_project_properties/sequencesPassingExclusion"
)

URL_EXPORT_REPLACEMENTS = (
    URL_PREFIX +
    "project/cov/module/covReplacementWebExporter"
)

URL_EXPORT_INSERTIONS = (
    URL_PREFIX +
    "project/cov/module/covInsertionWebExporter"
)

URL_EXPORT_DELETIONS = (
    URL_PREFIX +
    "project/cov/module/covDeletionWebExporter"
)

URL_DOWNLOAD = (
    URL_PREFIX +
    'glue_web_files/{webSubDirUuid}/{webFileName}'
)


@CACHE.memoize(expire=86400)
def invoke_download(url, function_name, gene, order, filename):
    click.echo('Request new {} for gene {}...'.format(filename, gene))
    resp = requests.post(
        url,
        json={
            'invoke-function': {
                'functionName': function_name,
                'argument': [
                    'CSV',
                    'variation.featureLoc.feature.name = {}'
                    .format(json.dumps(gene)),
                    order,
                    filename,
                    'LF'
                ]
            }
        }
    )
    resp.raise_for_status()
    payload = resp.json()
    click.echo('  Downloading {}...'.format(filename))
    url_download = URL_DOWNLOAD.format(**payload['tabularWebFileResult'])
    resp = requests.get(url_download)
    resp.raise_for_status()
    ret = resp.text
    click.echo('  Done.')
    return ret


@CACHE.memoize(expire=86400, typed=True)
def get_total_sequences():
    resp = requests.post(
        URL_NUM_SEQUENCES,
        json={'show': {'property': {'property': 'value'}}}
    )
    resp.raise_for_status()
    payload = resp.json()
    return int(payload['propertyValueResult']['value'])


def export_replacements(gene):
    return invoke_download(
        URL_EXPORT_REPLACEMENTS,
        'exportReplacements',
        gene,
        order=(
            '-num_seqs,+variation.featureLoc.feature.name,'
            '+codon_label_int,+replacement_aa'
        ),
        filename='replacements.csv'
    )


def export_insertions(gene):
    return invoke_download(
        URL_EXPORT_INSERTIONS,
        'exportInsertions',
        gene,
        order='-num_seqs,+last_ref_nt_before',
        filename='insertions.csv'
    )


def export_deletions(gene):
    return invoke_download(
        URL_EXPORT_DELETIONS,
        'exportDeletions',
        gene,
        order='-num_seqs,+reference_nt_start',
        filename='deletions.csv'
    )


def apply_napos_offset(napos, gene):
    for offset in GENE_OFFSETS[gene]:
        if napos >= offset['between'][0] and \
                napos <= offset['between'][1]:
            return napos - offset['offset']


def to_aapcnt(gene, replacements, insertions, deletions, total):
    rows = csv.DictReader(replacements.splitlines())
    aapcnts = []
    mutcounts = Counter()
    for row in rows:
        pos = int(row['codonNumber'])
        count = int(row['numSeqs'])
        mutcounts[pos] += count
        aapcnts.append({
            'gene': gene,
            'position': pos,
            'aa': row['replacementAminoAcid'],
            'percent': count / total,
            'count': count,
            'total': total
        })

    rows = csv.DictReader(deletions.splitlines())
    delcounts = Counter()
    for row in rows:
        # The GLUE data removes all NA deletions with length < 3,
        # however, they didn't fix the out-frame issue - which can
        # be addressed by codon alignment.
        count = int(row['numSeqs'])
        napos_start = int(row['refNtStartPosition'])
        napos_end = int(row['refNtEndPosition'])
        napos_start = apply_napos_offset(napos_start, gene)
        napos_end = apply_napos_offset(napos_end, gene)
        del_nalen = napos_end - napos_start + 1

        # deletion codon re-alignment start
        aapos_start = math.ceil((napos_start + 2) / 3)
        # deletion codon re-alignment end

        for pos in range(aapos_start, int(aapos_start + del_nalen / 3)):
            # del_nalen always >= 3
            mutcounts[pos] += count
            delcounts[pos] += count

    rows = csv.DictReader(insertions.splitlines())
    inscounts = Counter()
    for row in rows:
        # Similar to deletions, GLUE removes all NA insertions with
        # length < 3. Additional codon align adjustment is also
        # needed.
        pos = row['codonNumberBefore'].strip()
        if not pos:
            napos = int(row['refNtPositionBefore'])
            napos = apply_napos_offset(napos, gene)
            # insertion codon re-alignment start
            pos = math.floor((napos + 2) / 3)
            # insertion codon re-alignment end
        count = int(row['numSeqs'])
        pos = int(pos)
        mutcounts[pos] += count
        inscounts[pos] += count

    for pos0, aa in enumerate(GENEREFS[gene]):
        pos = pos0 + 1
        aapcnts.append({
            'gene': gene,
            'position': pos,
            'aa': aa,
            'percent': 1 - mutcounts[pos] / total,
            'count': total - mutcounts[pos],
            'total': total
        })
        if delcounts[pos] > 0:
            aapcnts.append({
                'gene': gene,
                'position': pos,
                'aa': 'd',
                'percent': delcounts[pos] / total,
                'count': delcounts[pos],
                'total': total
            })
        if inscounts[pos] > 0:
            aapcnts.append({
                'gene': gene,
                'position': pos,
                'aa': 'i',
                'percent': inscounts[pos] / total,
                'count': inscounts[pos],
                'total': total
            })
    return sorted(aapcnts, key=lambda p: (p['gene'], p['position'], p['aa']))


@click.command()
@click.argument('genes', nargs=-1)
def main(genes):
    totalseq = get_total_sequences()
    click.echo(
        'Total sequences passed exclusion criteria: {}'.format(totalseq)
    )
    os.makedirs(RESULTS_DIR, exist_ok=True)
    for gene in genes:
        reps = export_replacements(gene)
        ins = export_insertions(gene)
        dels = export_deletions(gene)
        destfile = os.path.join(RESULTS_DIR, 'aapcnt-{}.json'.format(gene))
        with open(destfile, 'w') as fp:
            payload = to_aapcnt(gene, reps, ins, dels, totalseq)
            json.dump(payload, fp, indent=2)
            print('Created {}'.format(destfile))


if __name__ == '__main__':
    main()
