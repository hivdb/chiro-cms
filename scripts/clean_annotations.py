#! /usr/bin/env python3

import re
import sys
import ruamel.yaml

yaml = ruamel.yaml.YAML()


def abort(msg, *args, **kwargs):
    print(msg.format(*args, **kwargs), file=sys.stderr)
    exit(1)


def valid_re(pattern):
    try:
        re.compile(pattern)
        return True
    except re.error:
        return False


def clean_annot_categories(data):
    categories = []
    dupnames = set()
    known_annots = {a['name'] for a in data['annotations']}

    for orig_cat in data['annotCategories']:
        name = orig_cat['name']
        if name in dupnames:
            abort('Duplicate annot category: {!r}', name)
        dupnames.add(name)
        display = orig_cat.get('display')
        dropdown = bool(orig_cat.get('dropdown', True))
        default_annot = orig_cat.get('defaultAnnot')
        default_annots = orig_cat.get('defaultAnnots', [])
        multi_select = bool(orig_cat.get('multiSelect'))
        annot_style = orig_cat.get('annotStyle')
        color = orig_cat.get('color')
        new_cat = {'name': name}
        if display is False or display:
            new_cat['display'] = display
        new_cat['dropdown'] = dropdown
        if multi_select:
            new_cat['defaultAnnots'] = [
                annot for annot in default_annots
                if annot in known_annots
            ]
        else:
            if default_annot is not None and \
                    default_annot not in known_annots:
                default_annot = None

            new_cat['defaultAnnot'] = default_annot
        new_cat['multiSelect'] = multi_select
        if annot_style:
            new_cat['annotStyle'] = annot_style
            if annot_style == 'aminoAcids' and color:
                new_cat['color'] = color
        else:
            new_cat['annotStyle'] = 'hide'
        categories.append(new_cat)
    data['annotCategories'] = categories


def clean_annotations(data):
    annotations = []
    dupnames = set()
    known_cats = {c['name'] for c in data['annotCategories']}

    for orig_annot in data['annotations']:
        name = orig_annot['name']
        label = orig_annot.get('label')
        level = orig_annot['level']
        category = orig_annot['category']
        if level == 'amino acid':
            level = 'aminoAcid'
        if name in dupnames:
            abort('Duplicated annot group found: {}', name)
        dupnames.add(name)
        if level not in ('position', 'aminoAcid'):
            abort('Invalid level {!r} for annot group {!r}', level, name)
        if category not in known_cats:
            abort('Unknown annot category: {!r}', category)
        hide_citations = bool(orig_annot.get('hideCitations'))
        color_rules = orig_annot.get('colorRules') or []
        if not all(valid_re(r) for r in color_rules):
            abort('Some color rules of annot group {!r} is invalid', name)
        new_annot = {'name': name}
        if label:
            new_annot['label'] = label
        new_annot.update({
            'level': level,
            'category': category,
            'hideCitations': hide_citations,
            'colorRules': color_rules
        })
        annotations.append(new_annot)
    data['annotations'] = annotations


def citation_id_sortkey(idx):
    id1, id2 = idx.split('.', 1)
    return (int(id1), int(id2))


def clean_citations(data):
    citations = []

    for cite in data['citations'].values():
        citations.append({
            'citationId': cite['citationId'],
            'sectionId': cite['sectionId'],
            'doi': cite['doi'],
            'author': cite['author'],
            'year': cite['year'],
            'section': cite['section']
        })
    citations.sort(key=lambda c: (c['citationId'], c['sectionId']))
    data['citations'] = {
        '{citationId}.{sectionId}'.format(**c): c
        for c in citations
    }


def clean_positions(data):
    positions = []
    known_annots = {a['name'] for a in data['annotations']}

    for posdata in data['positions']:
        position = posdata['position']
        annotations = posdata.get('annotations') or []
        old_aa_annots = posdata.get('aminoAcids') or []
        if old_aa_annots:
            new_aa_annots = {}
            for aadata in old_aa_annots:
                aa = aadata['aminoAcid']
                for annot in aadata['annotations']:
                    name = annot['name']
                    citation_ids = annot['citationIds']
                    new_aa_annots.setdefault(name, {
                        'name': name,
                        'aminoAcids': [],
                        'description': '',
                        'citationIds': set()
                    })
                    new_aa_annots[name]['aminoAcids'].append(aa)
                    new_aa_annots[name]['citationIds'] |= set(citation_ids)
            new_aa_annots = list(new_aa_annots.values())
            for annot in new_aa_annots:
                annot['citationIds'] = sorted(
                    annot['citationIds'],
                    key=citation_id_sortkey)
            annotations.extend(new_aa_annots)
            annotations.sort(key=lambda a: a['name'])
        annotations = [a for a in annotations if a['name'] in known_annots]
        if annotations:
            positions.append({
                'position': position,
                'annotations': annotations
            })
    data['positions'] = sorted(positions, key=lambda p: p['position'])


def main():
    if len(sys.argv) != 2:
        abort('Usage: {} <YAML_PATH>', sys.argv[0])
    yaml_path = sys.argv[1]
    with open(yaml_path) as read_fp:
        data = yaml.load(read_fp)
        clean_annot_categories(data)
        clean_annotations(data)
        clean_citations(data)
        clean_positions(data)
    data['version'] = '20200910120644'
    with open(yaml_path, 'w') as write_fp:
        yaml.dump(data, write_fp)


if __name__ == '__main__':
    main()