import os
import sys
import ruamel.yaml
from ruamel.yaml.comments import CommentedSeq

yaml = ruamel.yaml.YAML()


def main():
    if len(sys.argv) != 3:
        print('Usage: {} <INPUT_MUTANNOT_YAML> <OUTPUT_DIRECTORY>',
              file=sys.stderr)
        exit(1)

    file_in, dir_out = sys.argv[1:]
    os.makedirs(dir_out, exist_ok=True)
    with open(file_in) as fp:
        data = yaml.load(fp)
    for annotdef in data['annotations']:
        annot_name = annotdef['name']
        citations = {}
        positions = {}
        amino_acids = CommentedSeq()
        amino_acids.fa.set_flow_style()
        for posdata in data['positions']:
            pos = posdata['position']
            for annot in posdata['annotations']:
                if annot['name'] != annot_name:
                    continue
                for cid in annot['citationIds']:
                    cite = {**data['citations'][cid]}
                    cite.pop('citationId')
                    cite.pop('sectionId')
                    citations.setdefault(cid, {
                        **cite,
                        'positions': CommentedSeq()
                    })['positions'].append(pos)
                if annotdef['level'] == 'position':
                    positions.setdefault(annot['value'], {
                        'subgroup': annot['value'],
                        'description': annot.get('description', ''),
                        'positions': CommentedSeq()
                    })['positions'].append(pos)
                    if positions[annot['value']].get('description') == '':
                        positions[annot['value']].pop('description')
                else:
                    amino_acids.append(
                        '{}{}'.format(
                            pos, ''.join(sorted(annot['aminoAcids']))
                        )
                    )
        citations = list(citations.values())
        for cite in citations:
            cite['positions'].fa.set_flow_style()
        positions = list(positions.values())
        for subgroup in positions:
            subgroup['positions'].fa.set_flow_style()
        with open(os.path.join(
            dir_out, '{}.yml'.format(annot_name.lower().replace(' ', '-'))
        ), 'w') as fp:
            payload = {'name': annot_name}
            if 'label' in annotdef:
                payload['label'] = annotdef['label']
            payload.update({
                'level': annotdef['level'],
                'hideCitations': bool(annotdef.get('hideCitations')),
                'colorRules': annotdef.get('colorRules', []),
                'citations': citations
            })
            if annotdef['level'] == 'position':
                payload['positions'] = positions
            else:
                payload['aminoAcids'] = amino_acids
            yaml.dump(payload, fp)


if __name__ == '__main__':
    main()
