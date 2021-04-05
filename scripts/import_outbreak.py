import requests
import click
import ruamel.yaml
from pathlib import Path
import time
from datetime import datetime
from collections import defaultdict
from functools import lru_cache
import decimal
from decimal import Decimal
from operator import itemgetter

yaml = ruamel.yaml.YAML()

BASE_DIR = Path(__file__).absolute().parent.parent

RESULTS_DIR = BASE_DIR / 'resources' / 'outbreak-aapcnt'

API_MAIN = 'https://api.outbreak.info'
API_VERSION = API_MAIN + '/genomics/metadata'

API_ALL_MUTATIONS = API_MAIN + '/genomics/mutations?name=*'
API_GLOBAL_PREVALENCE_MUTATION = API_MAIN + \
        "/genomics/global-prevalence?&mutations={mutation}"

API_ALL_VARIANTS = API_MAIN + '/genomics/mutations-by-lineage'
API_GLOBAL_PREVALENCE_VARIANT = API_MAIN + \
        "/genomics/global-prevalence?" + \
        "pangolin_lineage={variant}"
# API_GLOBAL_PREVALENCE_VARIANT = API_MAIN + \
#         "/genomics/global-prevalence?" + \
#         "cumulative=true&pangolin_lineage={variant}"
API_VARIANT_MUTATIONS = API_MAIN + \
        "/genomics/lineage-mutations?pangolin_lineage={variant}&frequency=0.75"
# TODO: timestamp clean cache


def get_datetime_obj(datetime_str):
    return datetime.strptime(datetime_str, '%Y-%m-%d')


START_DATE = '2020-07-01'
START_DATE = get_datetime_obj(START_DATE)

FREQ_CALL_NUMBER = 5
FREQ_WAIT_SECONDS = 1
FREQ_CONTROL = 0


def freq_control(func):
    def wrapper(*args, **kwargs):
        global FREQ_CONTROL
        if FREQ_CONTROL >= FREQ_CALL_NUMBER:
            FREQ_CONTROL = 0
            time.sleep(FREQ_WAIT_SECONDS)
        FREQ_CONTROL += 1
        return func(*args, **kwargs)

    return wrapper


@freq_control
def query_api(api):
    resp = requests.get(api)
    resp = resp.json()
    return resp


def estimate_runtime(batch):
    total_count = len(batch)
    estimate_time = (total_count / FREQ_CALL_NUMBER) * (FREQ_WAIT_SECONDS + 1)

    print('Estimate run: {}s'.format(estimate_time))


decimal.getcontext().rounding = decimal.ROUND_HALF_UP


def round_number(float_number):
    if float_number > 1:
        return Decimal(str(float_number)).quantize(Decimal('1'))
    elif float_number > 0.1:
        return Decimal(str(float_number)).quantize(Decimal('1.0'))
    elif float_number < 0.01:
        return 0
    elif float_number == 0:
        return 0
    return Decimal(str(float_number)).quantize(Decimal('1.00'))


def get_proportion(numerator, denominator):
    proportion = round_number(
                numerator / denominator * 100)
    if proportion:
        proportion = '{}%'.format(proportion)
    else:
        return str(proportion)
    return proportion


def get_version():
    resp = query_api(API_VERSION)
    build_version = resp['build_version']
    date = build_version[:8]
    print('Build at: ' + date)


@lru_cache(maxsize=32)
def get_all_mutations():
    resp = query_api(API_ALL_MUTATIONS)

    mutations = resp['results']
    print('All mutations:', len(mutations))

    return mutations


def get_all_spike_mutations():
    all_mutations = get_all_mutations()

    spike_mutations = []
    for mut in all_mutations:
        name = mut['name']
        if name.startswith('s:'):
            spike_mutations.append(mut)

    print('All spike mutations:', len(spike_mutations))

    return spike_mutations


def get_mutation_prevalence():
    prevalence = defaultdict(list)

    all_spike_mutations = get_all_spike_mutations()

    estimate_runtime(all_spike_mutations)

    for mut in all_spike_mutations:
        name = mut['name']
        query = API_GLOBAL_PREVALENCE_MUTATION.format(mutation=name)

        resp = query_api(query)
        timepoints = resp['results']

        cummulative_lineage_count = 0
        cummulative_total_count = 0

        for tp in timepoints:
            date = tp['date']
            lineage_count = tp['lineage_count']
            cummulative_lineage_count += lineage_count
            total_count = tp['total_count']
            cummulative_total_count += total_count

            if get_datetime_obj(date) < START_DATE:
                continue

            if not date.endswith('01'):
                continue

            proportion = get_proportion(
                cummulative_lineage_count,
                cummulative_total_count
            )

            name = name.upper()
            prevalence[name].append({
                'date': date,
                'prevalence': proportion
            })

    prevalence = dict(prevalence)
    return prevalence


def process_mutations(save_dir):

    prevalence = get_mutation_prevalence()
    save_path = save_dir / 'aapcnt-mutations.yml'
    with open(save_path, 'w') as fp:
        yaml.dump(prevalence, fp)
        print('Updated {}'.format(save_path))


@lru_cache(maxsize=32)
def get_all_variants():
    resp = query_api(API_ALL_VARIANTS)

    variants = resp['results']
    print('All variants:', len(variants))

    variants_names = [i['pangolin_lineage'] for i in variants]

    return variants_names


def get_variant_mutations():
    mutations = {}

    all_variants = get_all_variants()
    estimate_runtime(all_variants)

    for variant in all_variants:
        variant = variant.upper()

        query = API_VARIANT_MUTATIONS.format(variant=variant)
        resp = query_api(query)

        mutations[variant] = []

        variant_mutations = defaultdict(list)
        for mut in resp['results']:
            gene = mut['gene']
            variant_mutations[gene].append({
                'gene': mut['gene'],
                'position': mut['codon_num'],
                'mutation': mut['mutation'].upper(),
            })

        variant_mutations = sorted(
            variant_mutations.items(),
            key=lambda x: x[0],
            reverse=True
        )

        for gene, muts in variant_mutations:
            muts.sort(key=itemgetter('position'))
            for mut in muts:
                mutations[variant].append(mut['mutation'])

        mutations[variant] = ','.join(mutations[variant])

    return mutations


def get_variant_global_prevalence():
    prevalence = defaultdict(list)

    all_variants = get_all_variants()
    estimate_runtime(all_variants)
    for variant in all_variants:
        query = API_GLOBAL_PREVALENCE_VARIANT.format(variant=variant)

        resp = query_api(query)
        timepoints = resp['results']

        cummulative_total_count = 0
        cummulative_lineage_count = 0

        for tp in timepoints:
            date = tp['date']
            lineage_count = tp['lineage_count']
            cummulative_lineage_count += lineage_count
            total_count = tp['total_count']
            cummulative_total_count += total_count

            if get_datetime_obj(date) < START_DATE:
                continue

            if not date.endswith('01'):
                continue

            proportion = get_proportion(
                cummulative_lineage_count,
                cummulative_total_count
            )

            variant = variant.upper()
            prevalence[variant].append({
                'date': date,
                'prevalence': proportion
            })
        # latest_tp = resp['results']

        # date = latest_tp['last_detected']
        # prevalence[variant].append({
        #     'date': date,
        #     'total_count': latest_tp['total_count'],
        #     'count': latest_tp['lineage_count'],
        #     'proportion': latest_tp['global_prevalence']
        # })

    prevalence = dict(prevalence)
    return prevalence


def process_variants(save_dir):

    prevalence = get_variant_global_prevalence()
    save_path = save_dir / 'aapcnt-variants.yml'
    with open(save_path, 'w') as fp:
        yaml.dump(prevalence, fp)
        print('Updated {}'.format(save_path))

    save_path = save_dir / 'variants-mutaions.yml'
    mutations = get_variant_mutations()
    with open(save_path, 'w') as fp:
        yaml.dump(mutations, fp)
        print('Updated {}'.format(save_path))


@click.command()
def import_outbreak():
    RESULTS_DIR.mkdir(exist_ok=True)

    get_version()

    # process_mutations(RESULTS_DIR)

    file_path = RESULTS_DIR / 'aapcnt-mutations.yml'
    with open(file_path) as fp:
        content = yaml.load(fp)


    # process_variants(RESULTS_DIR)


if __name__ == '__main__':
    import_outbreak()
