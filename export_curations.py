import os
import csv
import json
import pickle
import itertools
from collections import defaultdict
from indra.statements import *
from indra.statements import amino_acids
from indra.sources import indra_db_rest


header = ['ENTITYA', 'TYPEA', 'IDA', 'DATABASEA',
          'ENTITYB', 'TYPEB', 'IDB', 'DATABASEB',
          'EFFECT', 'MECHANISM', 'RESIDUE', 'SEQUENCE',
          'TAX_ID', 'CELL_DATA', 'TISSUE_DATA',
          'MODULATOR_COMPLEX', 'TARGET_COMPLEX',
          'MODIFICATIONA', 'MODASEQ', 'MODIFICATIONB', 'MODBSEQ',
          'PMID', 'DIRECT', 'NOTES', 'ANNOTATOR', 'SENTENCE']


def sanitize_text(text):
    if text is None:
        return text
    text = text.replace('[XREF_BIBR - XREF_BIBR]', '')
    text = text.replace('[XREF_BIBR]', '')
    text = text.replace('[XREF_BIBR', '')
    text = text.replace('XREF_BIBR', '')
    text = text.replace('(XREF_SUPPLEMENTARY)', '')
    text = text.strip()
    return text


comment_mapping = {
    'Ser15|EFFECT: down-regulates quantity by destabilization':
        'RESIDUE:S15;EFFECT:down-regulates quantity by destabilization',
    'Tyr15': 'RESIDUE:Y15',
    'Thr14': 'RESIDUE:T14',
    'Ser15': 'RESIDUE:S15',
    'effect:down-regulates activity': 'EFFECT:down-regulates activity',
    'effect:up-regulates activity': 'EFFECT:up-regulates activity',
    'residue:tyr705': 'RESIDUE:Y705',
    'direct:no': 'DIRECT:no',
    'TAXID:9606;TISSUE:BTO:0000759': 'TAXID:9606;CELL:BTO:0000759',
    'RESIDUE:S151;T753': 'RESIDUE:S151;RESIDUE:T753',
    'RESIDUE:Y448; SENTENCE:adaptor protein 3BP2 serves as a binding protein and a physiological substrate of SHP-1. 3BP2 is phosphorylated on tyrosyl residue 448 in response to TCR activation, and the phosphorylation is required for T c':
        'RESIDUE:Y448;SENTENCE:adaptor protein 3BP2 serves as a binding protein and a physiological substrate of SHP-1. 3BP2 is phosphorylated on tyrosyl residue 448 in response to TCR activation, and the phosphorylation is required for T c',
    'sentence:Integrin-bound PTP-PEST dephosphorylates RhoGDI1.':
        'SENTENCE:Integrin-bound PTP-PEST dephosphorylates RhoGDI1.',
    'EFFECT: dow-regulates quantity by degradation': 'EFFECT:down-regulates quantity by degradation',
}


def process_comment(comment):
    allowed_keys = {'CELL', 'TAXID', 'DIRECT', 'EFFECT', 'SENTENCE',
                    'MECHANISM', 'RESIDUE'}
    if not comment:
        return {}

    comment = comment_mapping.get(comment, comment)
    # Corner case: trailing ;
    parts = comment.strip(';').split(';')
    comment_data = defaultdict(list)
    for part in parts:
        part = part.strip()
        try:
            key, value = part.split(':', maxsplit=1)
        except ValueError:
            print('Value error when splitting comment part: %s' % part)
            print(comment)
            print('--')
            break
        if key not in allowed_keys:
            print('Key %s not in allowed keys' % key)
            print(comment)
            print('--')
            break
        if key == 'SENTENCE':
            comment_data['sentence'].append(value)
        elif key == 'EFFECT':
            comment_data['effect'].append(value)
        elif key == 'MECHANISM':
            comment_data['mechanism'].append(value)
        elif key == 'DIRECT':
            comment_data['direct'].append(value)
        elif key == 'CELL':
            comment_data['cell_data'].append(value)
        elif key == 'TAXID':
            comment_data['taxid'].append(value)
    return dict(comment_data)


def curations_to_rows(curations):
    # We need to check if there is activation and/or inhibition
    # among the statements. If there is, then we create a regulation
    # row.
    stmts_by_type = defaultdict(list)
    for stmt_package in curations:
        stmts_by_type[type(stmt_package[0])].append(stmt_package)
    has_pos_reg = (Activation in stmts_by_type
                   or IncreaseAmount in stmts_by_type)
    has_neg_reg = (Inhibition in stmts_by_type
                   or DecreaseAmount in stmts_by_type)
    has_mod = mod_class in stmts_by_type

    # Here we need to look at cases where there is a SENTENCE
    # involved

    # If there is no regulation statement we need to look at the
    # comment in the mod statement and see if there is a
    # mechanism that is a regulation statement. If there is, we
    # need to create a regulation statement with the mod statement
    # and the regulation statement.
    if not (has_pos_reg or has_neg_reg):
        if has_mod:
            for mod_stmt, mod_ev, mod_cur, mod_comment \
                    in stmts_by_type[mod_class]:
                if mod_comment and mod_comment.get('effect'):
                    # We just handle one effect for now
                    assert len(mod_comment['effect']) == 1
                    effect = mod_comment['effect'][0]
                    if 'down-regulates' in effect.lower():
                        stmt_cls = DecreaseAmount if \
                            effect.startswith('down-regulates quantity') else Inhibition
                        stmts_by_type[stmt_cls].append(
                            (stmt_cls(Agent('X'), Agent('Y')), None, None,
                             {'effect': [effect]})
                        )
                        has_neg_reg = True
                        mod_comment.pop('effect')
                    elif 'up-regulates' in effect.lower():
                        stmt_cls = IncreaseAmount if \
                            effect.startswith('up-regulates quantity') else Activation
                        stmts_by_type[stmt_cls].append(
                            (stmt_cls(Agent('X'), Agent('Y')), None, None,
                             {'effect': [effect]})
                        )
                        has_pos_reg = True
                        mod_comment.pop('effect')
                    else:
                        print('Unknown effect')
                        print(mod_stmt, mod_comment)
    # If there is no mechanism we need to look at the comment in the
    # regulation statement and see if there is a mechanism that is a mod
    # statement. If there is, we need to create a mod statement with the
    # regulation statement and the mod statement.
    elif not has_mod:
        for act_stmt, act_ev, act_cur, act_comment \
                in stmts_by_type[Activation] + stmts_by_type[Inhibition]:
            if act_comment and act_comment.get('mechanism') == [mod_key]:
                stmts_by_type[mod_class].append(
                    (mod_class(act_stmt.subj, act_stmt.obj), act_ev, act_cur, {})
                )
                has_mod = True
                act_comment.pop('mechanism')

    # If either mod or reg is missing that means that we couldn't
    # create any missing statements from comments and so have to
    # return
    if not has_mod or not (has_pos_reg or has_neg_reg):
        return []

    for mod_stmt_package, reg_stmt_package in \
            itertools.product(stmts_by_type[mod_class],
                              stmts_by_type[Activation] + stmts_by_type[Inhibition] +
                              stmts_by_type[IncreaseAmount] +
                              stmts_by_type[DecreaseAmount]):
        mod_stmt, mod_ev, mod_cur, mod_comment = mod_stmt_package
        reg_stmt, reg_ev, reg_cur, reg_comment = reg_stmt_package
        mod_comment_effect = mod_comment.get('effect')[0] if \
            mod_comment and mod_comment.get('effect') else None
        enz = mod_stmt.enz
        substrate = mod_stmt.sub
        is_activation = isinstance(reg_stmt, Activation)
        is_upergulation = isinstance(reg_stmt, IncreaseAmount)
        is_inhibition = isinstance(reg_stmt, Inhibition)
        is_downregulation = isinstance(reg_stmt, DecreaseAmount)

        if reg_comment and 'effect' in reg_comment:
            assert len(reg_comment['effect']) == 1
            effect = reg_comment['effect'][0].strip()
        elif mod_comment_effect:
            effect = mod_comment_effect.strip()
        else:
            if is_activation:
                effect = 'up-regulates'
            elif is_upergulation:
                effect = 'up-regulates quantity'
            elif is_inhibition:
                effect = 'down-regulates'
            elif is_downregulation:
                effect = 'down-regulates quantity'

        if mod_stmt.residue and mod_stmt.position:
            residue = \
                amino_acids[mod_stmt.residue]['short_name'].capitalize() + \
                mod_stmt.position
        else:
            residue = ''

        sentence_parts = []
        if mod_ev.text:
            sentence_parts.append(sanitize_text(mod_ev.text))
        if reg_ev and reg_ev.text:
            sentence_parts.append(sanitize_text(reg_ev.text))
        if mod_comment.get('sentence'):
            sentence_parts.extend(mod_comment['sentence'])
        if reg_comment and reg_comment.get('sentence'):
            sentence_parts.extend(reg_comment['sentence'])
        sentence = '|'.join(sorted(set(sentence_parts)))

        direct_comment = mod_comment.get('direct')
        if direct_comment and direct_comment[0].lower() == 'no':
            direct = 'NO'
        else:
            direct = 'YES'

        reg_taxid = reg_comment.get('taxid')
        mod_taxid = mod_comment.get('taxid')
        if reg_taxid:
            taxid = reg_taxid[0]
        elif mod_taxid:
            taxid = mod_taxid[0]
        else:
            taxid = '9606'

        curators = []
        if mod_cur and mod_cur.get('curator'):
            curators.append(mod_cur['curator'])
        if reg_cur and reg_cur.get('curator'):
            curators.append(reg_cur['curator'])
        parts = sorted(curators)[0].split('@')[0].split('.')
        curator = parts[0][0] + parts[1]

        yield [
            # 'ENTITYA', 'TYPEA', 'IDA', 'DATABASEA'
            enz.name, 'protein', enz.db_refs.get('UP'), 'UNIPROT',
            # 'ENTITYB', 'TYPEB', 'IDB', 'DATABASEB',
            substrate.name, 'protein', substrate.db_refs.get('UP'), 'UNIPROT',
            # 'EFFECT', 'MECHANISM', 'RESIDUE', 'SEQUENCE',
            effect, mod_key, residue, '',
            # 'TAX_ID', 'CELL_DATA', 'TISSUE_DATA',
            taxid, '', '',
            # 'MODULATOR_COMPLEX', 'TARGET_COMPLEX',
            '', '',
            # 'MODIFICATIONA', 'MODASEQ', 'MODIFICATIONB', 'MODBSEQ',
            '', '', '', '',
            # 'PMID', 'DIRECT', 'NOTES', 'ANNOTATOR'
            mod_ev.pmid, direct, '', curator,
            # 'SENTENCE'
            sentence
        ]


def merge_curation_rows(curation_rows):
    # Find rows without sites and find rows with sites
    # Merge sentences of rows without sites into rows with sites
    row_dicts = [dict(zip(header, row)) for row in curation_rows]
    rows_no_site = [row for row in row_dicts if not row['RESIDUE']]
    rows_with_site = [row for row in row_dicts if row['RESIDUE']]
    # We have to have both to do merging
    if rows_with_site and rows_no_site:
        new_rows = []
        for row in rows_no_site:
            for row_with_site in rows_with_site:
                if row['EFFECT'] == row_with_site['EFFECT']:
                    row_with_site_sentences = set(row_with_site['SENTENCE'].split('|'))
                    row_sentences = set(row['SENTENCE'].split('|'))
                    unique_sentences = '|'.join(list(row_sentences |
                                                     row_with_site_sentences))
                    new_row = row_with_site.copy()
                    new_row['SENTENCE'] = unique_sentences

                    new_rows.append(new_row)
    else:
        new_rows = row_dicts
    # Then look at the remaining rows and merge any that are exactly the same
    # except for the sentences, and merge the set of unique sentences into
    # one consolidated row
    rows_by_key = defaultdict(list)
    for row in new_rows:
        row_key = (row['EFFECT'], row['MECHANISM'], row['RESIDUE'])
        rows_by_key[row_key].append(row)
    merged_rows = []
    for key, rows in rows_by_key.items():
        if len(rows) == 1:
            merged_rows.append(rows[0].values())
        else:
            sentences = set()
            for row in rows:
                sentences |= set(row['SENTENCE'].split('|'))
            new_row = rows[0].copy()
            new_row['SENTENCE'] = '|'.join(list(sentences))
            merged_rows.append(new_row.values())
    return merged_rows


def get_pair_key(stmt, ev):
    subj = stmt.subj if isinstance(stmt, (RegulateActivity, RegulateAmount)) \
        else stmt.enz
    obj = stmt.obj if isinstance(stmt, (RegulateActivity, RegulateAmount)) \
        else stmt.sub
    return subj.name, obj.name, ev.pmid


if __name__ == '__main__':
    # -------------------------------
    #mod_key_short = 'dephos'
    #mod_key = 'dephosphorylation'
    #mod_class = Dephosphorylation
    mod_key_short = 'ubiq'
    mod_key = 'ubiquitination'
    mod_class = Ubiquitination
    # -------------------------------
    print('Exporting curations for %s' % mod_key)
    if os.path.exists('curations.json'):
        with open('curations.json', 'r') as fh:
            curs = json.load(fh)
    else:
        curs = indra_db_rest.get_curations()
        with open('curations.json', 'w') as fh:
            json.dump(curs, fh, indent=1)
    curs = [cur for cur in curs
            if cur.get('source') == 'signor_%s' % mod_key_short]
    # Sometimes we have duplicate curations that we can
    # squash here
    curs = {(cur['pa_hash'], cur['source_hash']): cur for cur in curs}.values()
    print('Found %d curations' % len(curs))
    with open('%ss_with_reg_sorted.pkl' % mod_key, 'rb') as fh:
        stmts = pickle.load(fh)
    stmts_by_hash = {stmt.get_hash(): stmt for stmt in stmts}
    ev_by_source_hash = {ev.source_hash: ev
                         for stmt in stmts for ev in stmt.evidence}
    # We need to get the statement corresponding to the curation
    # and the evidence corresponding to the curation
    # We also need to group curations so that the activity regulation
    # and the dephopshorylation are grouped and we can generate
    # a single curation row from it.

    stmts_by_pair_pubmed_key = defaultdict(list)

    curs_by_hashes = defaultdict(list)
    for cur in curs:
        if not cur['tag'] == 'correct':
            continue
        comment_data = process_comment(cur['text'])
        #if comment_data:
        #    print(comment_data)
        stmt = stmts_by_hash[cur['pa_hash']]
        ev = ev_by_source_hash[cur['source_hash']]
        curs_by_hashes[(stmt.get_hash(), ev.get_source_hash())].append(cur)
        key = get_pair_key(stmt, ev)
        stmts_by_pair_pubmed_key[key].append((stmt, ev, cur, comment_data))

    all_rows = []
    for key, curations in stmts_by_pair_pubmed_key.items():
        curation_rows = list(curations_to_rows(curations))
        curation_rows = merge_curation_rows(curation_rows)
        all_rows += curation_rows

    with open('%ss_with_reg_export.csv' % mod_key, 'wt') as fh:
        writer = csv.writer(fh, delimiter=',', quotechar='"')
        writer.writerow(header)
        for row in all_rows:
            writer.writerow(row)