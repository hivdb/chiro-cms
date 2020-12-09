import os
import csv
import json
import ruamel.yaml
from collections import defaultdict
from itertools import groupby, zip_longest, chain

from .func_mutannots import yield_mutannots_json

DRM_ANNOT_CATEGORY = 'escapeMutants'


def load_spike_mab_lookup(resource_dir):
    resource_name = os.path.join(resource_dir, 'mabs-table.yml')
    with open(resource_name, encoding='utf-8-sig') as fp:
        resdata = ruamel.yaml.load(fp, Loader=ruamel.yaml.Loader)
    lookup = defaultdict(list)
    for row in resdata:
        for ref in row['references']:
            refid = ref.get('refId')
            if refid:
                lookup[ref['refId']].extend(
                    ab['name'] for ab in row['antibodies']
                )
    return lookup


def build_citeid2refid_lookup(resource_dir, doi2citeid):
    """RefID is used by Bob

    This function create a lookup map convert CiteID to RefID
    """
    resource_name = os.path.join(resource_dir, 'refid_lookup.json')
    with open(resource_name, encoding='utf-8-sig') as fp:
        resdata = json.load(fp)
    lookup = {}
    for row in resdata:
        doi = row['doi']
        if doi in doi2citeid:
            lookup[doi2citeid[doi]] = row['refId']
    return lookup


def build_drm_citation_pairs(positions, drm_annot_names,
                             citeid2refid, citeid2doi,
                             all_citations):
    results = []
    for posdata in positions:
        pos = posdata['position']
        for annot in posdata['annotations']:
            if annot['name'] not in drm_annot_names:
                continue

            aas = annot['aminoAcids']
            for citeid in annot['citationIds']:
                citeid = int(citeid.split('.', 1)[0])
                if citeid not in citeid2refid:
                    raise KeyError(
                        'Unable to find RefID for citation {}. '
                        'Is it a preprint got published?'
                        .format(citeid2doi[citeid])
                    )
                results.append({
                    'position': pos,
                    'aminoAcids': ''.join(aas),
                    'refId': citeid2refid[citeid]
                })
    return results


def save_drm2citations(destpath, pairs):
    with open(destpath, 'w', encoding='utf-8-sig') as fp:
        writer = csv.DictWriter(fp, ['Refs', 'Pos', 'AAs'])
        writer.writeheader()
        pairs = sorted(pairs, key=lambda r: (r['position'], r['aminoAcids']))
        for (pos, aas), refids in groupby(
            pairs, lambda r: (r['position'], r['aminoAcids'])
        ):
            writer.writerow({
                'Pos': pos,
                'AAs': aas,
                'Refs': '; '.join(r['refId'] for r in refids)
            })
    print('create: {}'.format(destpath))


def save_citation2drms2mabs(destpath, pairs, refid2mabs):
    with open(destpath, 'w', encoding='utf-8-sig') as fp:
        writer = csv.DictWriter(fp, ['Ref', 'Mutations', 'MAbs'])
        writer.writeheader()
        pairs = sorted(pairs, key=lambda r: (r['refId']))
        for refid, muts in groupby(
            pairs, lambda r: (r['refId'])
        ):
            writer.writerow({
                'Ref': refid,
                'Mutations': '; '.join(
                    '{position}{aminoAcids}'.format(**mut)
                    for mut in muts
                ),
                'MAbs': ('\n'.join(refid2mabs[refid])
                         if refid2mabs[refid]
                         else 'NA')
            })
    print('create: {}'.format(destpath))


def save_mab2refs(destpath, refid2mabs):
    with open(destpath, 'w', encoding='utf-8-sig') as fp:
        writer = csv.DictWriter(fp, ['Refs', 'MAbs'])
        writer.writeheader()
        rows = [
            (refid, '\n'.join(mabs))
            for refid, mabs in refid2mabs.items()
        ]
        rows = sorted(rows, key=lambda r: r[1])
        writer.writerows({
            'MAbs': mabs,
            'Refs': '\n'.join(r[0] for r in refids)
        } for mabs, refids in groupby(rows, lambda r: r[1]))
    print('create: {}'.format(destpath))


def bobstyle_csvs(destpath, *srcpaths):
    all_rows = []
    num_cols = []
    for srcpath in srcpaths:
        with open(srcpath, encoding='utf-8-sig') as fp:
            rows = list(csv.reader(fp))
            all_rows.append(rows)
            num_cols.append(len(rows[0]) if rows else 0)
    with open(destpath, 'w', encoding='utf-8-sig') as fp:
        writer = csv.writer(fp)
        for row in zip_longest(*all_rows):
            row = list(chain(*[
                (r if r else [''] * num_col) + ['']
                for r, num_col in zip(row, num_cols)
            ]))
            if row:
                row.pop(-1)
            writer.writerow(row)
    print('create: {}'.format(destpath))


def build_drm_tables(resource_dir, build_dir, download_dir, **kw):
    for resname, payload, _ in yield_mutannots_json(resource_dir):
        all_citations = payload['citations']
        doi2citeid = {
            cite['doi']: cite['citationId']
            for cite in all_citations.values()
        }
        citeid2doi = {
            cite['citationId']: cite['doi']
            for cite in all_citations.values()
        }
        citeid2refid = build_citeid2refid_lookup(resource_dir, doi2citeid)

        if resname.lower() == 'spike':
            refid2mabs = load_spike_mab_lookup(resource_dir)
        else:
            refid2mabs = {}

        drm_annots = [
            annot for annot in payload['annotations']
            if annot['category'] == DRM_ANNOT_CATEGORY and
            annot['label'].lower() != 'all'
        ]
        drm_annot_names = {annot['name'] for annot in drm_annots}
        positions = [
            pos for pos in payload['positions']
            if any(
                annot['name'] in drm_annot_names
                for annot in pos['annotations']
            )
        ]
        drm_citation_pairs = build_drm_citation_pairs(
            positions, drm_annot_names,
            citeid2refid, citeid2doi,
            all_citations
        )
        save_drm2citations(
            os.path.join(download_dir, 'drms/{}-drm2refs.csv'.format(resname)),
            drm_citation_pairs
        )
        save_citation2drms2mabs(
            os.path.join(download_dir,
                         'drms/{}-ref2drms2mabs.csv'.format(resname)),
            drm_citation_pairs,
            refid2mabs
        )
        save_mab2refs(
            os.path.join(download_dir, 'drms/{}-mab2refs.csv'.format(resname)),
            refid2mabs
        )
        bobstyle_csvs(
            os.path.join(download_dir,
                         'drms/{}-merged-drms.csv'.format(resname)),
            os.path.join(download_dir,
                         'drms/{}-ref2drms2mabs.csv'.format(resname)),
            os.path.join(download_dir, 'drms/{}-drm2refs.csv'.format(resname)),
            os.path.join(download_dir, 'drms/{}-mab2refs.csv'.format(resname))
        )
