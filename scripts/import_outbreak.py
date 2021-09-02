#!/usr/bin/env python3
import requests
import click
import csv
import ruamel.yaml
from pathlib import Path
import simplejson
import time
from datetime import datetime
from collections import defaultdict
import re
import decimal
from decimal import Decimal
from operator import itemgetter

yaml = ruamel.yaml.YAML()

BASE_DIR = Path(__file__).absolute().parent.parent
RESULTS_DIR = BASE_DIR / 'resources' / 'outbreak-aapcnt'

FILE_CACHE_PATH = Path('/tmp') / 'py_file_cache'
FILE_CACHE_PATH.mkdir(exist_ok=True)


API_MAIN = 'https://api.outbreak.info'
API_VERSION = API_MAIN + '/genomics/metadata'

API_ALL_MUTATIONS = API_MAIN + '/genomics/mutations?name=*'
API_GLOBAL_PREVALENCE_MUTATION = API_MAIN + \
    "/genomics/global-prevalence?&mutations={mutation}"

API_ALL_VARIANTS = API_MAIN + '/genomics/lineage?name=*'
# API_ALL_VARIANTS = API_MAIN + '/genomics/mutations-by-lineage'
API_VARIANT_MUTATIONS = API_MAIN + \
    "/genomics/lineage-mutations?pangolin_lineage={variant}&frequency=0.75"
API_GLOBAL_PREVALENCE_VARIANT = API_MAIN + \
    "/genomics/global-prevalence?" + \
    "pangolin_lineage={variant}"
# API_GLOBAL_PREVALENCE_VARIANT = API_MAIN + \
#     "/genomics/global-prevalence?" + \
#     "cumulative=true&pangolin_lineage={variant}"

API_ALL_LOCATIONS = API_MAIN + \
    "/genomics/location?name=*"
API_VARIANT_GEOLOCATION = API_MAIN + \
    "/genomics/lineage-by-sub-admin-most-recent?" + \
    "pangolin_lineage={variant}"
# API_VARIANT_GEOLOCATION = API_MAIN + \
#     "/genomics/lineage-by-sub-admin-most-recent?" + \
#     "pangolin_lineage={variant}&detected=true"
API_VARIANT_USA_GEOLOCATION = API_MAIN + \
    "/genomics/lineage-by-sub-admin-most-recent?" + \
    "pangolin_lineage={variant}&detected=true&location_id=USA"
# API_VARIANT_LOCATION_PREV = API_MAIN + \
#     "/genomics/prevalence-by-location?" + \
#     "pangolin_lineage={variant}&location_id={location_id}&cumulative=true"


def get_datetime_obj(datetime_str):
    return datetime.strptime(datetime_str, '%Y-%m-%d')


START_DATE = get_datetime_obj('2020-01-01')

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
    try:
        resp = resp.json()
    except simplejson.errors.JSONDecodeError as e:
        print(api)
        print(resp.text)
        raise e
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


def map_with_skip(map_list, skip_list, operator):
    for rec in map_list:
        if rec in skip_list:
            continue
        result = operator(rec)

        yield rec, result


def file_cache(file_name):
    def wrapper(func):
        def wrapper2(*args, **kwargs):
            file_path = FILE_CACHE_PATH / file_name
            if file_path.exists():
                result = yaml.load(open(file_path))
            else:
                result = func(*args, **kwargs)
                with open(file_path, 'w') as fp:
                    yaml.dump(result, fp)

            return result

        return wrapper2

    return wrapper


def load_progress(progress_file_path):
    progress_path = FILE_CACHE_PATH / progress_file_path
    if not progress_path.exists():
        return None
    with open(progress_path) as fp:
        return yaml.load(fp)


def dump_progress(progress_list, progress_file_path):
    progress_path = FILE_CACHE_PATH / progress_file_path
    with open(progress_path, 'w') as fp:
        yaml.dump(progress_list, fp)


def dump_csv(file_path, records, headers=[]):
    if not records:
        return

    if headers:
        _records = []
        for r in records:
            _r = {}
            for h in headers:
                _r[h] = r.get(h, '')
            _records.append(_r)

        records = _records

    if not headers and records:
        headers = records[0].keys()

    for r in records:
        for k in r.keys():
            if k not in headers:
                del r[k]

    with open(file_path, 'w', encoding='utf-8-sig') as fd:
        writer = csv.DictWriter(fd, fieldnames=headers)
        writer.writeheader()
        writer.writerows(records)


def split_mut(mutation):
    if 'DEL' in mutation.upper():
        ref = ''
        start_pos, stop_pos = re.findall(r'(\d+)', mutation)
        pos = start_pos
        mut = 'del'
    else:
        ref, pos, mut = re.split(r'(\d+)', mutation)
        start_pos = pos
        stop_pos = ''
    return {
        'ref': ref,
        'pos': pos,
        'mut': mut,
        'start_pos': start_pos,
        'stop_pos': stop_pos,
    }


def get_proportion(numerator, denominator):
    proportion = round_number(
                numerator / denominator * 100)
    # if proportion:
    #     proportion = '{}%'.format(proportion)
    # else:
    #     return str(proportion)
    return str(proportion)


def filter_timepoints(timepoints):
    records = []
    MIN_PREVALENCE = 0.1
    for name, tp_list in timepoints.items():
        prev_list = [float(i['prevalence']) > MIN_PREVALENCE for i in tp_list]
        if not any(prev_list):
            continue
        rec = {
            'name': name
        }
        for tp in tp_list:
            date = tp['date']
            prevalence = tp['prevalence']
            rec[date] = prevalence
        records.append(rec)

    return records


def parse_mutation(mut_str):
    gene, mut_str = mut_str.split(':')
    pos_list = re.findall(r'\d+', mut_str)
    if len(pos_list) == 1:
        aa_start, aa_stop = pos_list[0], ''
    else:
        aa_start, aa_stop = pos_list

    if 'DEL' in mut_str.upper():
        ref = ''
        mut = 'del'
        name = '∆{}-{}'.format(aa_start, aa_stop)
    else:
        ref, _, mut = re.split(r'(\d+)', mut_str)
        name = '{}{}{}'.format(ref, aa_start, mut)

    return {
        'gene': gene,
        'ref': ref,
        'aa_start': aa_start,
        'aa_stop': aa_stop,
        'mut': mut,
        'ref_pos': '{}{}'.format(ref, aa_start),
        'name': name
    }


def get_version():
    resp = query_api(API_VERSION)
    build_version = resp['build_version']
    date = build_version[:8]
    print('Build at: ' + date)


@file_cache('mutations')
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
    all_spike_mutations = get_all_spike_mutations()
    all_spike_mutations = [m['name'] for m in all_spike_mutations]
    estimate_runtime(all_spike_mutations)

    processed_list = load_progress('processed_list_mut_tp_prev') or []
    prevalence = load_progress('mutation_tp_prev') or {}
    prevalence = defaultdict(list, prevalence)

    def operator(mutation):
        query = API_GLOBAL_PREVALENCE_MUTATION.format(mutation=mutation)
        resp = query_api(query)

        return resp

    for mutation, rec in \
            map_with_skip(all_spike_mutations, processed_list, operator):
        print(mutation)

        timepoints = rec['results']

        timepoints_group = defaultdict(list)

        for tp in timepoints:
            date = tp['date']
            year_month = date[:7]
            if get_datetime_obj(date) < START_DATE:
                continue

            timepoints_group[year_month].append(tp)

        prev_list = []

        for year_month, tp_list in timepoints_group.items():
            lineage_count = sum(tp['lineage_count'] for tp in tp_list)
            total_count = sum(tp['total_count'] for tp in tp_list)

            proportion = get_proportion(
                lineage_count,
                total_count
            )

            mutation_name = mutation.upper()
            prevalence[mutation_name].append({
                'date': year_month,
                'prevalence': proportion
            })

            prev_list.append(proportion)

        dump_progress(dict(prevalence), 'mutation_tp_prev')

        processed_list.append(mutation)
        dump_progress(processed_list, 'processed_list_mut_tp_prev')

        if all([float(p) < 0.1 for p in prev_list]):
            break

    prevalence = filter_timepoints(prevalence)
    return prevalence


@file_cache('all_variants')
def get_all_variants():
    resp = query_api(API_ALL_VARIANTS)

    variants = resp['results']
    print('All variants:', len(variants))

    variants_names = [i['name'] for i in variants]

    return variants_names


@file_cache('all_variants_count')
def get_all_variants_count():
    resp = query_api(API_ALL_VARIANTS)

    variants = {}
    for item in resp['results']:
        name = item['name'].upper()
        total_count = item['total_count']
        variants[name] = total_count
    return variants


def get_variant_mutations():
    all_variants = get_all_variants()
    variants = [i.upper() for i in all_variants]
    estimate_runtime(all_variants)

    processed_list = load_progress('processed_list') or []
    mutations = load_progress('variant_mutation_list') or {}

    def operator(variant):
        query = API_VARIANT_MUTATIONS.format(variant=variant)
        resp = query_api(query)
        print(query)

        return resp

    for variant, rec in map_with_skip(variants, processed_list, operator):
        print(variant)

        mutations[variant] = []

        variant_mutations = defaultdict(list)
        for name, mut_list in rec['results'].items():
            for mut in mut_list:
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

        mutations[variant] = ', '.join(mutations[variant])
        dump_progress(mutations, 'variant_mutation_list')

        processed_list.append(variant)
        dump_progress(processed_list, 'processed_list')

    records = []
    for variant, mut_list in mutations.items():
        records.append({
            'name': variant,
            'mutations': mut_list
        })

    return records


def get_variant_global_prevalence():
    all_variants = get_all_variants()
    variants = [i.upper() for i in all_variants]
    estimate_runtime(all_variants)

    processed_list = load_progress('processed_list_var_tp_prev') or []
    prevalence = load_progress('variant_tp_prev') or {}
    prevalence = defaultdict(list, prevalence)

    def operator(variant):
        query = API_GLOBAL_PREVALENCE_VARIANT.format(variant=variant)
        resp = query_api(query)

        return resp

    for variant, rec in map_with_skip(variants, processed_list, operator):

        timepoints = rec['results']

        variant = variant.upper()

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

            proportion = get_proportion(
                cummulative_lineage_count,
                cummulative_total_count
            )

            year_month = date[:7]

            prevalence[variant].append({
                'date': year_month,
                'prevalence': proportion
            })

        dump_progress(dict(prevalence), 'variant_tp_prev')

        processed_list.append(variant)
        dump_progress(processed_list, 'processed_list_var_tp_prev')

    prevalence = filter_timepoints(prevalence)
    return prevalence


@file_cache('locations')
def get_all_locations():
    resp = query_api(API_ALL_LOCATIONS)

    locations = resp['results']
    print('All locations:', len(locations))

    return locations


def get_variant_location_prevalence():
    all_variants = get_all_variants()
    variants = [i.upper() for i in all_variants]
    estimate_runtime(all_variants)

    processed_list = load_progress('processed_list_var_loc_prev') or []
    prevalence = load_progress('variant_loc_prev') or {}
    prevalence = defaultdict(list, prevalence)

    def operator(variant):
        query = API_VARIANT_GEOLOCATION.format(variant=variant)
        resp = query_api(query)

        return resp

    for variant, rec in map_with_skip(variants, processed_list, operator):
        variant = variant.upper()
        print(variant)

        results = rec['results']
        for loc in results:
            proportion = loc['proportion'] * 100
            proportion = str(round_number(proportion))

            location = loc['name']
            prevalence[variant].append({
                'location': location,
                'prevalence': proportion
            })

        dump_progress(dict(prevalence), 'variant_loc_prev')

        processed_list.append(variant)
        dump_progress(processed_list, 'processed_list_var_loc_prev')

    return dict(prevalence)


def collect_variant_mutations(save_dir):
    save_path = save_dir / 'variants-mutations.yml'
    mutations = get_variant_mutations()
    with open(save_path, 'w') as fp:
        yaml.dump(mutations, fp)
        print('Updated {}'.format(save_path))

    records = []
    for item in mutations:
        name = item['name']
        mut_list = [
            i.strip() for i in item['mutations'].split(',') if i.strip()]
        for mut in mut_list:
            gene, mut = mut.split(':')
            mut = split_mut(mut)
            rec = {
                'name': name,
                'gene': gene,
            }
            rec.update(mut)
            records.append(rec)
    dump_csv(save_dir / 'variants-mutations.csv', records)


def make_variant_mut_tracking_table(save_dir, spike_only=True):
    save_path = save_dir / 'variants-mut-tracking.csv'
    variants_mut = get_variant_mutations()

    variants_count = get_all_variants_count()

    for rec in variants_mut:
        name = rec['name']
        rec['total_count'] = variants_count[name]

    s_mutations = defaultdict(int)
    other_mutations = defaultdict(int)
    for rec in variants_mut:
        s_mut_ident = []
        for mut in rec['mutations'].split(','):
            mut = mut.upper().strip()
            if not mut:
                continue

            mut = parse_mutation(mut)
            if mut['gene'] == 'S':
                s_mut_ident.append(mut)
                s_mutations[mut['ref_pos']] += 1
                rec[mut['ref_pos']] = mut['mut']
            else:
                other_mutations[mut['ref_pos']] += 1
                rec[mut['ref_pos']] = mut['mut']
        del rec['mutations']

        s_mut_ident.sort(key=itemgetter('aa_start', 'ref'))

        rec['spike_ident'] = ','.join([s['name'] for s in s_mut_ident])

    s_mut = sorted([{
        'name': k,
        'count': v,
        }
        for k, v in s_mutations.items()], key=lambda x: -x['count'])
    s_mut = [s['name'] for s in s_mut]

    o_mut = sorted([{
        'name': k,
        'count': v,
        }
        for k, v in other_mutations.items()], key=lambda x: -x['count'])
    o_mut = [s['name'] for s in o_mut]

    if not spike_only:
        all_mut = s_mut + o_mut
        records = variants_mut
    else:
        all_mut = s_mut
        spike_ident_groups = defaultdict(list)
        for i in variants_mut:
            spike_ident = i['spike_ident']
            spike_ident_groups[spike_ident].append(i)

        records = []
        for spike_ident, var_list in spike_ident_groups.items():
            pangolins = ', '.join(sorted([v['name'] for v in var_list]))
            r = {
                'ident': spike_ident,
                'pangolins': pangolins,
                'total_count': sum(v['total_count'] for v in var_list)
            }
            for k, v in var_list[0].items():
                if k == 'total_count':
                    continue
                r[k] = v

            records.append(r)

    header = ['ident', 'pangolins'] + ['total_count'] + all_mut

    dump_csv(save_path, records, header)
    print('Updated {}'.format(save_path))


def collect_variant_time_prevalence(save_dir):
    prevalence = get_variant_global_prevalence()
    save_path = save_dir / 'aapcnt-variants.yml'
    with open(save_path, 'w') as fp:
        yaml.dump(prevalence, fp)
        print('Updated {}'.format(save_path))

    dump_csv(
        RESULTS_DIR / 'aapcnt-variants.csv',
        prevalence
    )


def collect_variant_location_prevalence(save_dir):
    prevalence = get_variant_location_prevalence()
    save_path = save_dir / 'aapcnt-loc-variants.yml'
    with open(save_path, 'w') as fp:
        yaml.dump(prevalence, fp)
        print('Updated {}'.format(save_path))


def collect_mutation_time_prevalence(save_dir):
    prevalence = get_mutation_prevalence()
    save_path = save_dir / 'aapcnt-mutations.yml'
    with open(save_path, 'w') as fp:
        yaml.dump(prevalence, fp)
        print('Updated {}'.format(save_path))

    headers = list(prevalence[0].keys())
    headers.insert(0, 'Mut')
    headers.insert(0, 'Pos')
    headers.insert(0, 'Ref')
    headers.insert(0, 'Gene')

    for rec in prevalence:
        mutation = rec['name']
        gene, mutation = mutation.split(':', 1)

        if 'DEL' in mutation.upper():
            ref = ''
            pos = re.search(r'(\d+)', mutation).group()
            mut = 'del'
        else:
            ref, pos, mut = re.split(r'(\d+)', mutation)

        rec['Gene'] = gene
        rec['Ref'] = ref
        rec['Pos'] = pos
        rec['Mut'] = mut
    dump_csv(
        RESULTS_DIR / 'aapcnt-mutations.csv',
        prevalence,
        headers,
    )


@click.command()
def import_outbreak():
    RESULTS_DIR.mkdir(exist_ok=True)

    get_version()

    collect_variant_mutations(RESULTS_DIR)
    # make_variant_mut_tracking_table(RESULTS_DIR)

    # collect_variant_time_prevalence(RESULTS_DIR)

    # collect_variant_location_prevalence(RESULTS_DIR)

    # collect_mutation_time_prevalence(RESULTS_DIR)


if __name__ == '__main__':
    import_outbreak()
