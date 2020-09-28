import os
import re
import json
import ruamel.yaml

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


def build_mutannots_json(resource_dir, buildres_dir, **kw):
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
                annotdata_lookup[annotdata['name']] = annotdata
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
                    cite_pos = cite.pop('positions')
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
                        for pos in subgroup['positions']:
                            posdata = positions.setdefault(pos, {
                                'position': pos,
                                'annotations': []
                            })
                            posdata['annotations'].append({
                                'name': annot_name,
                                'value': subgroup['subgroup'],
                                'description': subgroup.get('description', ''),
                                'citationIds': pos_cites[pos]
                            })
                else:
                    for aas in annotdata['aminoAcids']:
                        match = re.match(r'^(\d+)([A-Zid*]+)$', aas)
                        pos, aas = match.groups()
                        pos = int(pos)
                        aas = list(aas)
                        posdata = positions.setdefault(pos, {
                            'position': pos,
                            'annotations': []
                        })
                        posdata['annotations'].append({
                            'name': annot_name,
                            'description': '',
                            'aminoAcids': aas,
                            'citationIds': pos_cites[pos]
                        })

        payload = {
            'taxonomy': geneconfig['taxonomy'],
            'gene': geneconfig['gene'],
            'refSequence': geneconfig['refSequence'],
            'fragmentOptions': geneconfig['fragmentOptions'],
            'annotCategories': categories,
            'annotations': annotdefs,
            'citations': citations,
            'positions': list(positions.values()),
            'version': VERSION
        }
        dest_json = os.path.join(
            buildres_dir, 'mutannot-{}.json'.format(resname))
        with open(dest_json, 'w') as fp:
            json.dump(payload, fp, indent=2)
            print('create: {}'.format(dest_json))
