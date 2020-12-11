import os
import csv
import json
from itertools import groupby, zip_longest, chain

from .func_mutannots import yield_mutannots_json

DRM_ANNOT_CATEGORY = 'escapeMutants'


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


def build_triplets(
    positions, drm_annot_names,
    all_citations, resource_dir
):
    results = []
    citeid2refid = build_citeid2refid_lookup(
        resource_dir, all_citations
    )
    citeid2doi = {
        cite['citationId']: cite['doi']
        for cite in all_citations.values()
    }
    for posdata in positions:
        pos = posdata['position']
        for annot in posdata['annotations']:
            if annot['name'] not in drm_annot_names:
                continue

            aas = annot['aminoAcids']
            aa_attrs = annot.get('aminoAcidAttrs', {})
            for citeid in annot['citationIds']:
                citeid = int(citeid.split('.', 1)[0])
                if citeid not in citeid2refid:
                    raise KeyError(
                        'Unable to find RefID for citation {}. '
                        'Perhaps a preprint just got published?'
                        .format(citeid2doi[citeid])
                    )
                for aa in aas:
                    for mab in aa_attrs.get(aa, {}).get('resistance', []):
                        results.append({
                            'position': pos,
                            'aminoAcid': aa,
                            'refId': citeid2refid[citeid],
                            'mAb': mab
                        })
    return results


def sort_groupby(items, key):
    items = sorted(items, key=key)
    return groupby(items, key)


def uniq_join_attrs(items, attrname, by=''):
    return by.join(sorted({item[attrname] for item in items}))


def save_triplets(destpath, triplets):
    with open(destpath, 'w', encoding='utf-8-sig') as fp:
        writer = csv.DictWriter(fp, ['position', 'aminoAcid', 'refId', 'mAb'])
        writer.writeheader()
        for triplet in triplets:
            triplet['mAb'] = '="{mAb}"'.format(**triplet)
        writer.writerows(triplets)
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


def build_drm_tables(resource_dir, build_dir, download_dir, **kw):
    for resname, payload, _ in yield_mutannots_json(resource_dir):
        all_citations = payload['citations']

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
        triplets = build_triplets(
            positions, drm_annot_names,
            all_citations, resource_dir
        )

        dest_triplets = os.path.join(
            download_dir, 'drms/{}-triplets.csv'.format(resname)
        )
        save_triplets(dest_triplets, triplets)

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
