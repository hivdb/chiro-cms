import os
import re
import warnings
import ruamel.yaml  # type: ignore

yaml = ruamel.yaml.YAML()

VERSION = '20200924115632'


def get_citation_id(citation, reverse_citations):
    doi = citation.get('doi')
    refid = citation.get('refID')
    section = citation['section']
    reverse_citation = reverse_citations.setdefault(doi or refid, {
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
    ranges = obj.get('positions', [])
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
    if not positions:
        positions = (pos for pos, _ in get_amino_acids(obj))
    return sorted(set(positions))


def get_amino_acids(obj):
    results = set()
    all_posaas = obj.get('aminoAcids', [])
    for posaas in all_posaas:
        pos, aas = parse_mut(posaas)
        for aa in aas:
            results.add((pos, aa))
    return sorted(results)


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
        raise ValueError(
            'Invalid mutation string: {!r}; use "i" or "d" for indels'
            .format(mut)
        )
    pos, aas = match.groups()
    pos = int(pos)
    if '_' in aas:
        aas = [aas]
    else:
        aas = list(aas)
    return pos, aas


def build_aa_attrs(annotdata, citationId, sectionId):
    all_posaas = get_amino_acids(annotdata)
    attrs = annotdata.get('aminoAcidAttrs', [])
    lookup = {}
    for attrdef in attrs:
        attr = attrdef['attr']
        value = attrdef.get('value')
        pairs = attrdef.get('pairs', [])
        if isinstance(pairs, dict):
            pairs = list(pairs.items())
        for mut, attrval in pairs:
            pos, aas = parse_mut(mut)
            for aa in aas:
                if (pos, aa) not in lookup:
                    lookup[(pos, aa)] = {}
                if attr in lookup[(pos, aa)]:
                    warnings.warn(
                        'Unexpected behavior may be observed when duplicated '
                        'attribute {!r} was defined for {}{}'
                        .format(attr, pos, aa)
                    )
                lookup[(pos, aa)][attr] = attrval
        if value:
            for pos, aa in all_posaas:
                if (pos, aa) not in lookup:
                    lookup[(pos, aa)] = {}
                if attr not in lookup[(pos, aa)]:
                    lookup[(pos, aa)][attr] = value
                else:
                    warnings.warn(
                        'Unexpected behavior may be observed when "value" and '
                        '"pairs" are both defined for attribute {!r}'
                        .format(attr)
                    )
    grouped_by_pos = {}
    for (pos, aa), attrs in lookup.items():
        grouped_by_pos.setdefault(pos, []).append({
            'aminoAcid': aa,
            'citationId': citationId,
            'sectionId': sectionId,
            **attrs
        })
    return grouped_by_pos


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
                all_posaas = []
                all_posaa_attrs = {}

                for cite in annotdata['citations']:
                    cite_pos = get_positions(cite)
                    cite_id = get_citation_id(cite, rev_citations)
                    cite_id_str = '{citationId}.{sectionId}'.format(**cite_id)
                    citations[cite_id_str] = {
                        **cite_id,
                        'author': cite['author'],
                        'year': cite['year'],
                        'doi': cite.get('doi'),
                        'refID': cite.get('refID'),
                        'section': cite['section']
                    }
                    for pos in cite_pos:
                        pos_cites.setdefault(pos, []).append(cite_id_str)
                    if level == 'aminoAcid':
                        all_posaas.extend(get_amino_acids(cite))
                        for pos, attrs in build_aa_attrs(cite,
                                                         **cite_id).items():
                            all_posaa_attrs.setdefault(pos, []).extend(attrs)
                all_posaas = sorted(set(all_posaas))

                if level == 'position':
                    if 'positions' not in annotdata:
                        raise KeyError(
                            "'positions' is required for annotation {!r}"
                            .format(annot_name)
                        )
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
                    for pos, aa in all_posaas:
                        posdata = positions.setdefault(pos, {
                            'position': pos,
                            'annotations': {}
                        })
                        annot = posdata['annotations'].setdefault(annot_name, {
                            'name': annot_name,
                            'description': '',
                            'aminoAcids': [],
                            'aminoAcidAttrs': [],
                            'citationIds': pos_cites[pos]
                        })
                        annot['aminoAcids'].append(aa)
                    for pos in set(pos for pos, _ in all_posaas):
                        annot = positions[pos]['annotations'][annot_name]
                        annot['aminoAcidAttrs'].extend(
                            all_posaa_attrs.get(pos, [])
                        )
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
