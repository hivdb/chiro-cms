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
    keywords = func_def.get('keywords', {})
    func = func_mapping[func]
    return func, args, keywords


def cond_not(subcond_def, posdata, payload):
    func, args, kw = get_func_args(subcond_def, CONDITION_FUNCS)
    return not func(*args, **kw, posdata=posdata, payload=payload)


def cond_all(*subcond_defs, posdata, payload):
    for subcond_def in subcond_defs:
        func, args, kw = get_func_args(subcond_def, CONDITION_FUNCS)
        flag = func(*args, **kw, posdata=posdata, payload=payload)
        if not flag:
            return False
    return True


def readable_aas(aas):
    return [aa.replace('d', 'del').replace('i', 'ins')
            for aa in aas]


def cond_any(*subcond_defs, posdata, payload):
    for subcond_def in subcond_defs:
        func, args, kw = get_func_args(subcond_def, CONDITION_FUNCS)
        flag = func(*args, **kw, posdata=posdata, payload=payload)
        if flag:
            return True
    return False


def cond_has_annot(annot_name, *, posdata, payload, subgroup=None):
    for annot in posdata['annotations']:
        if annot['name'] == annot_name:
            if subgroup is None or annot['value'] == subgroup:
                return True
    return False


def cond_has_no_annot(annot_name, *, posdata, payload, subgroup=None):
    for annot in posdata['annotations']:
        if annot['name'] == annot_name:
            if subgroup is None or annot['value'] == subgroup:
                return False
    return True


def cond_has_annot_except_subgroup(
    annot_name, *subgroups, posdata, payload
):
    for annot in posdata['annotations']:
        if annot['name'] == annot_name and annot['value'] not in subgroups:
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
                aas='/'.join(readable_aas(sorted(aas))),
                be='are' if is_plural else 'is'
            ),
            is_plural=is_plural
        )
    raise ValueError(
        'Position {} does not have annotation {}'.format(pos, annot_name)
    )


def prefix_mutation_of_annotcat(cat_name, *, posdata, payload):
    annot_names = set()
    for annot in payload['annotations']:
        if annot['category'] != cat_name:
            continue
        annot_names.add(annot['name'])

    refseq = payload['refSequence']
    pos = posdata['position']
    aas = set()
    for annot in posdata['annotations']:
        if annot['name'] not in annot_names:
            continue
        aas |= set(annot['aminoAcids'])
    if not aas:
        raise ValueError(
            'Field aminoAcids of position {}, category {} is empty',
            pos, cat_name
        )
    is_plural = len(aas) > 1
    return Prefix(
        text='{refaa}{pos}{aas} {be}'.format(
            refaa=refseq[pos - 1],
            pos=pos,
            aas='/'.join(readable_aas(sorted(aas))),
            be='are' if is_plural else 'is'
        ),
        is_plural=is_plural
    )


def triggered_aas_mutation(annot_name, *, posdata, payload):
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
        return ''.join(aas)
    raise ValueError(
        'Position {} does not have annotation {}'.format(pos, annot_name)
    )


def triggered_aas_mutation_of_annot_cat(cat_name, *, posdata, payload):
    pos = posdata['position']
    annot_names = set()
    for annot in payload['annotations']:
        if annot['category'] != cat_name:
            continue
        annot_names.add(annot['name'])
    aas = set()
    for annot in posdata['annotations']:
        if annot['name'] not in annot_names:
            continue
        aas |= set(annot['aminoAcids'])
    if not aas:
        raise ValueError(
            'Field aminoAcids of position {}, category {} is empty',
            pos, cat_name
        )
    return ''.join(aas)


def variable_get_citations(annot_name, *, posdata, payload):
    citations = []
    all_citations = payload['citations']
    for annot in posdata['annotations']:
        if annot['name'] != annot_name:
            continue
        citation_ids = annot['citationIds']
        for cid in citation_ids:
            cite = all_citations[cid]
            if cite.get('doi'):
                citations.append('[^{doi}]'.format(**cite))
            else:
                citations.append('[^{refID}]'.format(**cite))
        break
    return ''.join(citations)


def variable_get_subgroup(annot_name, *, posdata, payload):
    pos = posdata['position']
    for annot in posdata['annotations']:
        if annot['name'] != annot_name:
            continue
        return annot['value']
    raise ValueError(
        'Position {} does not have annotation {}'.format(pos, annot_name)
    )


LAST_PUNC_WILDCARD = '\x13\x13LAST_PUNC_WILDCARD\x13\x13'


def natlang_join(items, footnotes=None, *,
                 conj='and',
                 last_punc=LAST_PUNC_WILDCARD, oxford_comma=False):
    size = len(items)
    if size == 0:
        return ''
    if not footnotes:
        footnotes = [''] * size
    buffers = []
    if LAST_PUNC_WILDCARD in items[0]:
        buffers.append(items[0])
    else:
        buffers.append('{}{}'.format(items[0], LAST_PUNC_WILDCARD))
    if size > 1:
        for item, fn in zip(items[1:-1], footnotes[:-2]):
            buffers[-1] = buffers[-1].replace(LAST_PUNC_WILDCARD, ',')
            if LAST_PUNC_WILDCARD in item:
                buffers.append('{} {}'.format(fn, item))
            else:
                buffers.append('{} {}{}'.format(fn, item, LAST_PUNC_WILDCARD))
        if oxford_comma:
            buffers[-1] = buffers[-1].replace(LAST_PUNC_WILDCARD, ',')
        else:
            buffers[-1] = buffers[-1].replace(LAST_PUNC_WILDCARD, '')
        if LAST_PUNC_WILDCARD in items[-1]:
            buffers.append('{} {} {}'.format(footnotes[-2], conj, items[-1]))
        else:
            buffers.append('{} {} {}{}'.format(
                footnotes[-2], conj, items[-1], LAST_PUNC_WILDCARD
            ))
    buffers[-1] = buffers[-1].replace(LAST_PUNC_WILDCARD, last_punc)
    buffers.append(footnotes[-1])
    return ''.join(buffers)


def get_prefix_by_gram_num(length, singular, plural):
    if length <= 1:
        prefix = singular
    else:
        prefix = plural
    if prefix:
        prefix = prefix + ' '
    else:
        prefix = ''
    return prefix


def variable_join_annot_cat_subgroup(
    cat_name, *, conj='and', footnote=False,
    singular_prefix=None,
    plural_prefix=None,
    posdata, payload
):
    pos = posdata['position']
    subgroups = []
    footnotes = []
    for annot in payload['annotations']:
        if annot['category'] != cat_name:
            continue
        for posannot in posdata['annotations']:
            if posannot['name'] == annot['name']:
                subgroups.append(posannot['value'])
                if footnote:
                    footnotes.append(variable_get_citations(
                        annot['name'], posdata=posdata, payload=payload
                    ))
                else:
                    footnotes.append('')
                break
    if len(subgroups) > 0:
        prefix = get_prefix_by_gram_num(
            len(subgroups), singular_prefix, plural_prefix)

        return prefix + natlang_join(subgroups, footnotes, conj=conj)
    else:
        raise ValueError(
            'Position {} does not have any annotation belong to {}'
            .format(pos, cat_name)
        )


def variable_join_annot_label_of_annot_cat(
    cat_name, *, conj='and', footnote=False,
    singular_prefix=None,
    plural_prefix=None,
    posdata, payload
):
    pos = posdata['position']
    labels = []
    footnotes = []
    for annot in payload['annotations']:
        if annot['category'] != cat_name:
            continue
        for posannot in posdata['annotations']:
            if posannot['name'] == annot['name']:
                labels.append(annot['label'])
                if footnote:
                    footnotes.append(variable_get_citations(
                        annot['name'], posdata=posdata, payload=payload
                    ))
                else:
                    footnotes.append('')
                break
    if len(labels) > 0:
        prefix = get_prefix_by_gram_num(
            len(labels), singular_prefix, plural_prefix)
        return prefix + natlang_join(labels, footnotes, conj=conj)
    else:
        raise ValueError(
            'Position {} does not have any annotation belong to {}'
            .format(pos, cat_name)
        )


def variable_join_aa_attr(
    *, annot_name, attr, format='{}', conj='and', posdata, payload
):
    pos = posdata['position']
    result = []
    for annot in posdata['annotations']:
        if annot['name'] != annot_name:
            continue
        for one in annot['aminoAcidAttrs']:
            value = one.get(attr)
            result.append(format.format(value))
    if not result:
        raise ValueError(
            'Position {} does not have annotation {}'.format(pos, annot_name)
        )
    return natlang_join(result, conj=conj)


CONDITION_FUNCS = {
    'all': cond_all,
    'any': cond_any,
    'not': cond_not,
    'hasAnnot': cond_has_annot,
    'hasAnnotExceptSubgroup': cond_has_annot_except_subgroup,
    'hasAnnotExceptSubgroups': cond_has_annot_except_subgroup,
    'hasAnnotCategory': cond_has_annot_category
}


PREFIX_FUNCS = {
    'position': prefix_position,
    'mutation': prefix_mutation,
    'mutationOfAnnotCat': prefix_mutation_of_annotcat
}


VARIABLE_FUNCS = {
    'getCitations': variable_get_citations,
    'getSubgroup': variable_get_subgroup,
    'joinSubgroupOfAnnotCat': variable_join_annot_cat_subgroup,
    'joinAnnotLabelOfAnnotCat': variable_join_annot_label_of_annot_cat,
    'joinAminoAcidAttr': variable_join_aa_attr
}

TRIGGERED_AAS_FUNCS = {
    'mutation': triggered_aas_mutation,
    'mutationOfAnnotCat': triggered_aas_mutation_of_annot_cat
}


def check_condition(cond_def, posdata, payload):
    func, args, kw = get_func_args(cond_def, CONDITION_FUNCS)
    return func(*args, **kw, posdata=posdata, payload=payload)


def format_prefix(prefix_def, posdata, payload):
    func, args, kw = get_func_args(prefix_def, PREFIX_FUNCS)
    return func(*args, **kw, posdata=posdata, payload=payload)


def extract_variables(variable_defs, posdata, payload):
    variables = {}
    for key, variable_def in variable_defs.items():
        func, args, kw = get_func_args(variable_def, VARIABLE_FUNCS)
        variables[key] = func(*args, **kw, posdata=posdata, payload=payload)
    return variables


def extract_triggered_aas(triggered_aas_def, posdata, payload):
    if not triggered_aas_def:
        return set()
    if triggered_aas_def == '*':
        return set('ACDEFGHIKLMNPQRSTVWYid')
    func, args, kw = get_func_args(triggered_aas_def, TRIGGERED_AAS_FUNCS)
    return set(
        func(*args, **kw, posdata=posdata, payload=payload)
    )


def build_references_markdown(payload):
    buffers = []
    for citation in payload['citations'].values():
        if citation.get('doi'):
            buffers.append(
                '[^{doi}]: {author} {year}, '
                '[doi.org/{doi}](https://doi.org/{doi}).'
                .format(**citation)
            )
        else:
            buffers.append(
                '[^{refID}]: {author} {year}.'
                .format(**citation)
            )

    return '\n'.join(buffers)


def build_mutannots_comments(resource_dir, buildres_dir, **kw):
    for resname, payload, geneconfig in yield_mutannots_json(resource_dir):
        cond_comments = sorted(geneconfig.get('conditionalComments', []),
                               key=lambda cc: cc['rank'])

        all_comments = []
        refseq = payload['refSequence']

        for posdata in payload['positions']:
            should_include = False
            position = posdata['position']
            refaa = refseq[position - 1]
            par_comments = defaultdict(list)
            par_footnotes = defaultdict(list)
            triggered_aas = set()

            for cmt_def in cond_comments:
                needother = cmt_def.get('requireOtherComment', False)
                cond_def = cmt_def['condition']
                if not check_condition(cond_def, posdata, payload):
                    continue
                if not needother:
                    should_include = True
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
                variables = cmt_def.get('variables', {})
                variables = extract_variables(variables, posdata, payload)
                tpl = Template(tpl)
                text = tpl.substitute(variables)
                par_comments[prefix.text].append(text)

                fn_tpl = cmt_def.get('footnote')
                if fn_tpl:
                    fn_tpl = Template(fn_tpl)
                    fn_text = fn_tpl.substitute(variables)
                    par_footnotes[prefix.text].append(fn_text)
                else:
                    par_footnotes[prefix.text].append('')
                triggered_aas |= extract_triggered_aas(
                    cmt_def.get('triggeredAAs'), posdata, payload)

            pos_comments = []
            for prefix, all_text in par_comments.items():
                all_fn = par_footnotes[prefix]
                comment_text = natlang_join(
                    all_text, all_fn, conj='and',
                    last_punc='.', oxford_comma=True)
                pos_comments.append('{} {}'.format(prefix, comment_text))

            if should_include and pos_comments:
                if refaa in triggered_aas:
                    triggered_aas.remove(refaa)
                all_comments.append({
                    'position': position,
                    'triggeredAAs': ''.join(sorted(triggered_aas)),
                    'comment': ' '.join(pos_comments)
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
                },  {
                    # 'name': 'triggeredAAs',
                    # 'label': 'AAs'
                    # }, {
                    'name': 'comment',
                    'label': 'Comment',
                    'textAlign': 'left',
                    'sortable': False
                }],
                'data': all_comments,
                'references': build_references_markdown(payload)
            }, fp, indent=2)
            print('create: {}'.format(dest_json))
