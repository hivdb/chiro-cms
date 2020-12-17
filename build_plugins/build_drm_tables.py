import os
import csv
import json
from itertools import groupby, zip_longest, chain

from .func_mutannots import yield_mutannots_json
from .generefs import get_refaa

EPITOPE_ANNOT_CATEGORY = 'epitopes'
DRM_ANNOT_CATEGORY = 'resistance'


def build_citeid2refid_lookup(resource_dir, all_citations):
    """RefID is used by Bob

    This function create a lookup map convert CiteID to RefID
    """
    doi2citeid = {
        cite['doi']: cite['citationId']
        for cite in all_citations.values()
    }
    resource_name = os.path.join(resource_dir, 'refid_lookup.json')
    with open(resource_name, encoding='utf-8-sig') as fp:
        resdata = json.load(fp)
    lookup = {}
    for row in resdata:
        doi = row['doi']
        if doi in doi2citeid:
            lookup[doi2citeid[doi]] = row['refId']
    return lookup


def build_epitope_lookup(positions, epitope_annot_names):
    lookup = {}
    suffix_len = len(' epitope')
    for posdata in positions:
        pos = posdata['position']
        pos_lookup = lookup[pos] = set()
        for annot in posdata['annotations']:
            annot_name = annot['name']
            if annot_name not in epitope_annot_names:
                continue
            mab = annot_name[:-suffix_len]
            pos_lookup.add(mab)
    return lookup


def get_triplets(
    positions, drm_annot_names,
    all_citations, resource_dir,
    epitope_lookup, mabs_with_epitope
):
    results = []
    citeid2refid = build_citeid2refid_lookup(
        resource_dir, all_citations
    )
    citeid2doi = {
        cite['citationId']: cite['doi']
        for cite in all_citations.values()
    }
    annot_prefix_len = len('resistance ')
    for posdata in positions:
        pos = posdata['position']
        for annot in posdata['annotations']:
            if annot['name'] not in drm_annot_names:
                continue
            mab = annot['name'][annot_prefix_len:]
            all_attrs = annot['aminoAcidAttrs']
            for attrs in all_attrs:
                citeid = attrs['citationId']
                if citeid not in citeid2refid:
                    raise KeyError(
                        'Unable to find RefID for citation {}. '
                        'Perhaps a preprint just got published?'
                        .format(citeid2doi[citeid])
                    )
                results.append({
                    **attrs,
                    'position': pos,
                    'refId': citeid2refid[attrs['citationId']],
                    'mAb': mab,
                    'hasEpitopeData': mab in mabs_with_epitope,
                    'isEpitopeAny': bool(epitope_lookup[pos]),
                    'isEpitopeSelf': mab in epitope_lookup[pos]
                })

    return sorted(
        results,
        key=lambda t: (
            t['position'],
            t['aminoAcid'],
            t['refId'],
            t['mAb']
        )
    )


def sort_groupby(items, key):
    items = sorted(items, key=key)
    return groupby(items, key)


def uniq_join_attrs(items, attrname, by=''):
    return by.join(sorted({item[attrname] for item in items}))


def save_triplets(destpath, triplets, gene):
    with open(destpath, 'w', encoding='utf-8-sig') as fp:
        writer = csv.DictWriter(
            fp, ['position', 'refAA', 'aminoAcid', 'refId',
                 'mAb', 'studyType', 'log10Fold', 'hasEpitopeData',
                 'isEpitopeSelf', 'isEpitopeAny'],
            extrasaction='ignore'
        )
        writer.writeheader()
        for triplet in triplets:
            pos = triplet['position']
            writer.writerow({
                **triplet,
                'refAA': get_refaa(gene, pos),
                'mAb': '={}'.format(json.dumps(triplet['mAb'])),
                'hasEpitopeData': 'Yes' if triplet['hasEpitopeData'] else '',
                'isEpitopeSelf': 'Yes' if triplet['isEpitopeSelf'] else '',
                'isEpitopeAny': 'Yes' if triplet['isEpitopeAny'] else ''
            })
    print('create: {}'.format(destpath))


def save_drm2citations(destpath, triplets):
    with open(destpath, 'w', encoding='utf-8-sig') as fp:
        writer = csv.DictWriter(fp, ['Pos', 'AAs', 'Refs'])
        writer.writeheader()
        for pos, partials in sort_groupby(triplets, lambda r: r['position']):
            partials = list(partials)

            aas = uniq_join_attrs(partials, 'aminoAcid')
            refids = uniq_join_attrs(partials, 'refId', '; ')

            writer.writerow({
                'Pos': pos,
                'AAs': aas,
                'Refs': refids
            })
    print('create: {}'.format(destpath))


def save_ref2drms2mabs(destpath, triplets):
    with open(destpath, 'w', encoding='utf-8-sig') as fp:
        writer = csv.DictWriter(fp, ['Ref', 'Mutations', 'MAbs'])
        writer.writeheader()
        for refid, partials in sort_groupby(triplets, lambda r: r['refId']):
            partials = list(partials)

            muts = '; '.join(
                '{}{}'.format(pos, uniq_join_attrs(aas, 'aminoAcid'))
                for pos, aas in sort_groupby(partials, lambda t: t['position'])
            )
            mabs = uniq_join_attrs(partials, 'mAb', '; ')

            writer.writerow({
                'Ref': refid,
                'Mutations': muts,
                'MAbs': mabs if mabs else 'NA'
            })
    print('create: {}'.format(destpath))


def save_mabs2refs(destpath, triplets):
    with open(destpath, 'w', encoding='utf-8-sig') as fp:
        writer = csv.DictWriter(fp, ['Refs', 'MAbs'])
        writer.writeheader()
        rows = [
            (refid, uniq_join_attrs(mabs, 'mAb', '; '))
            for refid, mabs in
            sort_groupby(triplets, lambda r: r['refId'])
        ]
        rows = sort_groupby(rows, lambda r: r[1])
        writer.writerows({
            'MAbs': mabs,
            'Refs': '; '.join(r[0] for r in refids)
        } for mabs, refids in rows)
    print('create: {}'.format(destpath))


def save_drm2mabs(destpath, triplets):
    with open(destpath, 'w', encoding='utf-8-sig') as fp:
        writer = csv.DictWriter(fp, ['Pos', 'AAs', 'MAbs'])
        writer.writeheader()
        for pos, partials in sort_groupby(triplets, lambda r: r['position']):
            partials = list(partials)

            aas = uniq_join_attrs(partials, 'aminoAcid')
            mabs = uniq_join_attrs(partials, 'mAb', '; ')

            writer.writerow({
                'Pos': pos,
                'AAs': aas,
                'MAbs': mabs
            })
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


def get_category_annots(payload, category):
    return [
        annot for annot in payload['annotations']
        if annot['category'] == category
    ]


def build_drm_tables(resource_dir, build_dir, download_dir, **kw):
    for resname, payload, _ in yield_mutannots_json(resource_dir):
        gene = 'S'
        if resname == 'rdrp':
            gene = 'RdRP'
        all_citations = payload['citations']

        drm_annots = get_category_annots(payload, DRM_ANNOT_CATEGORY)
        epitope_annots = get_category_annots(payload, EPITOPE_ANNOT_CATEGORY)
        drm_annot_names = {annot['name'] for annot in drm_annots}
        epitope_annot_names = {annot['name'] for annot in epitope_annots}

        positions = [
            pos for pos in payload['positions']
            if any(
                annot['name'] in drm_annot_names
                for annot in pos['annotations']
            )
        ]
        epitope_lookup = build_epitope_lookup(
            positions, epitope_annot_names)
        mabs_with_epitope = set(chain(*epitope_lookup.values()))
        triplets = get_triplets(
            positions, drm_annot_names,
            all_citations, resource_dir,
            epitope_lookup, mabs_with_epitope
        )

        dest_triplets = os.path.join(
            download_dir, 'drms/{}-triplets.csv'.format(resname)
        )
        save_triplets(dest_triplets, triplets, gene)

        dest_drm2refs = os.path.join(
            download_dir, 'drms/{}-drm2refs.csv'.format(resname)
        )
        save_drm2citations(dest_drm2refs, triplets)

        dest_ref2drms2mabs = os.path.join(
            download_dir, 'drms/{}-ref2drms2mabs.csv'.format(resname)
        )
        save_ref2drms2mabs(dest_ref2drms2mabs, triplets)

        dest_mabs2refs = os.path.join(
            download_dir, 'drms/{}-mabs2refs.csv'.format(resname)
        )
        save_mabs2refs(dest_mabs2refs, triplets)

        dest_drm2mabs = os.path.join(
            download_dir, 'drms/{}-drm2mabs.csv'.format(resname)
        )
        save_drm2mabs(dest_drm2mabs, triplets)

        bobstyle_csvs(
            os.path.join(download_dir,
                         'drms/{}-merged-drms.csv'.format(resname)),
            dest_ref2drms2mabs,
            dest_drm2refs,
            dest_drm2mabs,
            dest_mabs2refs
        )
