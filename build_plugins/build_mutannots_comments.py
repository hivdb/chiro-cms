import os
import json

from string import Template
from collections import defaultdict, namedtuple

from .func_mutannots import yield_mutannots_json


Prefix = namedtuple('Prefix', ['text', 'is_plural'])
PartialComment = namedtuple('PartialComment', ['prefix', 'text'])


def get_func_args(func_def, func_mapping):
    func = func_def['func']
    args = func_def.get('args', [])
    func = func_mapping[func]
    return func, args


def cond_has_annot(annot_name, *, posdata, payload):
    for annot in posdata['annotations']:
        if annot['name'] == annot_name:
            return True
    return False


def cond_has_annot_category(cat_name, *, posdata, payload):
    for annot in payload['annotations']:
        if annot['category'] != cat_name:
            continue
        if cond_has_annot(annot['name'], posdata=posdata, payload=payload):
            return True
    return False


def prefix_position(*, posdata, payload):
    refseq = payload['refSequence']
    pos = posdata['position']
    return Prefix(
        text='{refaa}{pos} is'.format(
            refaa=refseq[pos - 1], pos=pos),
        is_plural=False
    )


def prefix_mutation(annot_name, *, posdata, payload):
    refseq = payload['refSequence']
    pos = posdata['position']
    for annot in posdata['annotations']:
        if annot['name'] != annot_name:
            continue
        aas = annot['aminoAcids']
        if not aas:
            raise ValueError(
                'Field aminoAcids of position {}, annotation {} is empty',
                pos, annot_name
            )
        is_plural = len(aas) > 1
        return Prefix(
            text='{refaa}{pos}{aas} {be}'.format(
                refaa=refseq[pos - 1],
                pos=pos,
                aas='/'.join(sorted(aas)),
                be='are' if is_plural else 'is'
            ),
            is_plural=is_plural
        )
    raise ValueError(
        'Position {} does not have annotation {}'.format(pos, annot_name)
    )


def value_get_subgroup(annot_name, *, posdata, payload):
    pos = posdata['position']
    for annot in posdata['annotations']:
        if annot['name'] != annot_name:
            continue
        return annot['value']
    raise ValueError(
        'Position {} does not have annotation {}'.format(pos, annot_name)
    )


def natlang_join(items, *, oxford_comma=False):
    if len(items) > 1:
        return (
            ', '.join(items[:-1]) +
            (',' if oxford_comma else '') +
            ' and ' +
            items[-1]
        )
    elif len(items) == 1:
        return items[0]
    else:
        return ''


def value_join_annot_cat_subgroup(cat_name, *, posdata, payload):
    pos = posdata['position']
    subgroups = []
    for annot in payload['annotations']:
        if annot['category'] != cat_name:
            continue
        for posannot in posdata['annotations']:
            if posannot['name'] == annot['name']:
                subgroups.append(posannot['value'])
                break
    if len(subgroups) > 0:
        return natlang_join(subgroups)
    else:
        raise ValueError(
            'Position {} does not have any annotation belong to {}'
            .format(pos, cat_name)
        )


CONDITION_FUNCS = {
    'hasAnnot': cond_has_annot,
    'hasAnnotCategory': cond_has_annot_category
}


PREFIX_FUNCS = {
    'position': prefix_position,
    'mutation': prefix_mutation
}


VALUE_FUNCS = {
    'getSubgroup': value_get_subgroup,
    'joinAnnotCatSubgroup': value_join_annot_cat_subgroup
}


def check_condition(cond_def, posdata, payload):
    func, args = get_func_args(cond_def, CONDITION_FUNCS)
    return func(*args, posdata=posdata, payload=payload)


def format_prefix(prefix_def, posdata, payload):
    func, args = get_func_args(prefix_def, PREFIX_FUNCS)
    return func(*args, posdata=posdata, payload=payload)


def extract_values(value_defs, posdata, payload):
    values = {}
    for key, value_def in value_defs.items():
        func, args = get_func_args(value_def, VALUE_FUNCS)
        values[key] = func(*args, posdata=posdata, payload=payload)
    return values


def build_mutannots_comments(resource_dir, buildres_dir, **kw):
    for resname, payload, geneconfig in yield_mutannots_json(resource_dir):
        cond_comments = sorted(geneconfig.get('conditionalComments', []),
                               key=lambda cc: cc['rank'])

        all_comments = []

        for posdata in payload['positions']:
            position = posdata['position']
            par_comments = defaultdict(list)

            for cmt_def in cond_comments:
                cond_def = cmt_def['condition']
                if not check_condition(cond_def, posdata, payload):
                    continue
                prefix = cmt_def['prefix']
                prefix = format_prefix(prefix, posdata, payload)

                if prefix.is_plural:
                    tpl = cmt_def.get(
                        'templatePlural', cmt_def.get('template')
                    )
                else:
                    tpl = cmt_def.get(
                        'templateSingular', cmt_def.get('template')
                    )
                if not tpl:
                    raise ValueError(
                        'Template of {commentName} is empty'.format(**cmt_def)
                    )
                values = cmt_def.get('values', {})
                values = extract_values(values, posdata, payload)
                tpl = Template(tpl)
                text = tpl.substitute(values)
                par_comments[prefix.text].append(text)

            pos_comments = []
            for prefix, all_text in par_comments.items():
                comment_text = natlang_join(all_text, oxford_comma=True)
                pos_comments.append('{} {}'.format(prefix, comment_text))

            if pos_comments:
                all_comments.append({
                    'position': position,
                    'comment': '. '.join(pos_comments) + '.'
                })
        all_comments = sorted(all_comments, key=lambda c: c['position'])
        dest_json = os.path.join(
            buildres_dir, '{}-comments.json'.format(resname))
        with open(dest_json, 'w') as fp:
            json.dump({
                'columnDefs': [{
                    'name': 'position',
                    'label': 'Pos',
                    'sort': 'numeric'
                }, {
                    'name': 'comment',
                    'label': 'Comment',
                    'textAlign': 'left',
                    'sortable': False
                }],
                'data': all_comments
            }, fp, indent=2)
            print('create: {}'.format(dest_json))
