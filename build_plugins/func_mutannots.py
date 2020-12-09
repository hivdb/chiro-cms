import os
import re
import warnings
import ruamel.yaml
from collections import defaultdict

yaml = ruamel.yaml.YAML()

VERSION = '20200924115632'


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


def load_config_and_data(resource_dir):
    mutannots_dir = os.path.join(resource_dir, 'mutannots')
    for resyaml in os.listdir(mutannots_dir):
        resyaml_path = os.path.join(mutannots_dir, resyaml)
        match = re.search(r'(.+)-config\.ya?ml$', resyaml)
        if not os.path.isfile(resyaml_path) or not match:
            continue
        with open(resyaml_path, encoding='utf-8-sig') as fp:
            geneconfig = yaml.load(fp)

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
                if isinstance(annotdata, dict):
                    annotdata = [annotdata]
                for ad in annotdata:
                    annotdata_lookup[ad['name']] = ad

        yield resname, geneconfig, annotdata_lookup


def parse_mut(mut):
    match = re.match(
        r'^[A-Z]?(\d+)([A-Zid*]+|[A-Z*]_[A-Z*]+)$', mut)
    if not match:
        raise ValueError('Invalid mutation string: {!r}'.format(mut))
    pos, aas = match.groups()
    pos = int(pos)
    if '_' in aas:
        aas = [aas]
    else:
        aas = list(aas)
    return pos, aas


def build_aa_attrs_lookup(annotdata):
    attrs = annotdata.get('aminoAcidAttrs', [])
    lookup = defaultdict(dict)
    for attrdef in attrs:
        attr = attrdef['attr']
        pairs = attrdef['pairs']
        if isinstance(pairs, dict):
            pairs = list(pairs.items())
        for mut, attrval in pairs:
            pos, aas = parse_mut(mut)
            for aa in aas:
                if attr in lookup[(pos, aa)]:
                    warnings.warn(
                        'Duplicated attribute {!r} for {}{}'
                        .format(attr, pos, aa)
                    )
                lookup[(pos, aa)][attr] = attrval
    return lookup


def yield_mutannots_json(resource_dir):
    for resname, geneconfig, annotdata_lookup in \
            load_config_and_data(resource_dir):
        categories = []
        annotdefs = []
        citations = {}
        positions = {}
        rev_citations = {}
        for cat in geneconfig['annotCategories']:
            cat_annots = cat.pop('annotNames')
            categories.append(cat)
            for annot_name in cat_annots:
                pos_cites = {}
                annotdata = annotdata_lookup[annot_name]
                level = annotdata['level']
                annotdef = {
                    'name': annot_name,
                    'label': annotdata.get('label'),
                    'level': level,
                    'category': cat['name'],
                    'hideCitations': bool(annotdata.get('hideCitations')),
                    'colorRules': annotdata.get('colorRules', [])
                }
                if not annotdef['label']:
                    annotdef.pop('label')
                annotdefs.append(annotdef)

                for cite in annotdata['citations']:
                    cite_pos = get_positions(cite)
                    cite_id = get_citation_id(cite, rev_citations)
                    cite_id_str = '{citationId}.{sectionId}'.format(**cite_id)
                    citations[cite_id_str] = {
                        **cite_id,
                        'author': cite['author'],
                        'year': cite['year'],
                        'doi': cite['doi'],
                        'section': cite['section']
                    }
                    for pos in cite_pos:
                        pos_cites.setdefault(pos, []).append(cite_id_str)

                if level == 'position':
                    for subgroup in annotdata['positions']:
                        for pos in get_positions(subgroup):
                            posdata = positions.setdefault(pos, {
                                'position': pos,
                                'annotations': {}
                            })
                            posdata['annotations'][annot_name] = {
                                'name': annot_name,
                                'value': subgroup['subgroup'],
                                'description': subgroup.get('description', ''),
                                'citationIds': pos_cites[pos]
                            }
                else:
                    attr_lookup = build_aa_attrs_lookup(annotdata)
                    for aas in annotdata['aminoAcids']:
                        pos, aas = parse_mut(aas)
                        posdata = positions.setdefault(pos, {
                            'position': pos,
                            'annotations': {}
                        })
                        annot = posdata['annotations'].setdefault(annot_name, {
                            'name': annot_name,
                            'description': '',
                            'aminoAcids': [],
                            'aminoAcidAttrs': {},
                            'citationIds': pos_cites[pos]
                        })
                        annot['aminoAcids'].extend(aas)
                        annot['aminoAcidAttrs'].update({
                            aa: attr_lookup[(pos, aa)]
                            for aa in aas
                            if (pos, aa) in attr_lookup
                        })
        positions = list(positions.values())
        for posdata in positions:
            posdata['annotations'] = list(posdata['annotations'].values())

        payload = {
            'taxonomy': geneconfig['taxonomy'],
            'gene': geneconfig['gene'],
            'refSequence': geneconfig['refSequence'],
            'fragmentOptions': geneconfig['fragmentOptions'],
            'annotCategories': categories,
            'annotations': annotdefs,
            'citations': citations,
            'positions': positions,
            'version': VERSION
        }

        yield resname, payload, geneconfig
