import os
import re
import json
import ruamel.yaml

yaml = ruamel.yaml.YAML()

VERSION = '20200910120644'


def get_citation_id(citation, reverse_citations):
    doi = citation['doi']
    section = citation['section']
    reverse_citation = reverse_citations.setdefault(doi, {
        'idx': len(reverse_citations) + 1,
        'sections': {}
    })
    citation_id = reverse_citation['idx']
    section_lookup = reverse_citation['sections']
    section_id = section_lookup.setdefault(section, 1)
    return {
        'citationId': citation_id,
        'sectionId': section_id
    }


def get_positions(obj):
    positions = []
    ranges = obj.get('positions')
    for rangetext in ranges:
        if isinstance(rangetext, int):
            positions.append(rangetext)
        elif rangetext.isdigit():
            positions.append(int(rangetext))
        else:
            start, end = rangetext.split('-', 1)
            start = int(start.strip())
            end = int(end.strip())
            positions.extend(range(start, end + 1))
    return sorted(set(positions))


def build_mutannots_table(resource_dir, buildres_dir, **kw):
    mutannots_dir = os.path.join(resource_dir, 'mutannots')
    for resyaml in os.listdir(mutannots_dir):
        resyaml_path = os.path.join(mutannots_dir, resyaml)
        match = re.search(r'(.+)-config\.ya?ml$', resyaml)
        if not os.path.isfile(resyaml_path) or not match:
            continue
        with open(resyaml_path, encoding='utf-8-sig') as fp:
            geneconfig = yaml.load(fp)
        refseq = geneconfig['refSequence']

        resname = match.group(1)
        genedir = os.path.join(mutannots_dir, resname)
        annotdata_lookup = {}
        for annotyaml in os.listdir(genedir):
            match = re.search(r'.+\.ya?ml$', annotyaml)
            if not match:
                continue
            with open(os.path.join(genedir, annotyaml),
                      encoding='utf-8-sig') as fp:
                annotdata = yaml.load(fp)
                annotdata_lookup[annotdata['name']] = annotdata
        coldefs = [
            {'name': 'position',
             'label': 'Pos',
             'sort': 'numeric'},
            {'name': 'refAA',
             'label': 'Ref'}
        ]
        rows = {
            pos0 + 1: {
                'position': pos0 + 1,
                'refAA': aa
            } for pos0, aa in enumerate(refseq)
        }
        for cat in geneconfig['annotCategories']:
            cat_annots = cat.pop('annotNames')
            for annot_name in cat_annots:
                annotdata = annotdata_lookup[annot_name]
                annot_label = annotdata.get(
                    'tableColLabel',
                    annotdata.get('label', annot_name)
                )
                level = annotdata['level']

                pos_by_cite = [
                    set(get_positions(cite)) for cite in annotdata['citations']
                ]

                if all(pos == pos_by_cite[0] for pos in pos_by_cite):
                    pos_by_cite = pos_by_cite[:1]
                    col_names = [annot_label]
                    coldefs.append({
                        'name': annot_label,
                        'label': annot_label,
                        'none': ''
                    })
                else:
                    col_names = [
                        '{}-{author}-{year}'.format(annot_label, **cite)
                        for cite in annotdata['citations']
                    ]
                    coldefs.extend([{
                        'name': ('{}-{author}-{year}'
                                 .format(annot_label, **cite)),
                        'label': ('{}\n({}{})'
                                  .format(annot_label,
                                          cite['author'],
                                          cite['year'] - 2000)),
                        'none': ''
                    } for cite in annotdata['citations']])

                for col_name, posset in zip(col_names, pos_by_cite):
                    if level == 'position':
                        for subgroup in annotdata['positions']:
                            for pos in get_positions(subgroup):
                                if pos not in posset:
                                    continue
                                row = rows[pos]
                                row[col_name] = subgroup['subgroup']
                    else:
                        for aas in annotdata['aminoAcids']:
                            match = re.match(
                                r'^(\d+)([A-Zid*]+|[A-Z*]_[A-Z*]+)$', aas)
                            pos, aas = match.groups()
                            pos = int(pos)
                            if pos not in posset:
                                continue
                            row = rows[pos]
                            if '_' in aas:
                                row[col_name] = aas
                            else:
                                row[col_name] = ''.join(sorted(aas))

        rows = list(rows.values())
        dest_json = os.path.join(
            buildres_dir, 'mutannot-{}-table.json'.format(resname))
        with open(dest_json, 'w') as fp:
            json.dump({
                'columnDefs': coldefs,
                'data': rows
            }, fp, indent=2)
            print('create: {}'.format(dest_json))
