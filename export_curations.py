import os
import csv
import json
import pickle
import itertools
from collections import defaultdict
from indra.statements import Dephosphorylation, RegulateActivity, \
    Activation, Inhibition, Agent
from indra.statements import amino_acids
from indra.sources import indra_db_rest


header = ['ENTITYA', 'TYPEA', 'IDA', 'DATABASEA',
          'ENTITYB', 'TYPEB', 'IDB', 'DATABASEB',
          'EFFECT', 'MECHANISM', 'RESIDUE', 'SEQUENCE',
          'TAX_ID', 'CELL_DATA', 'TISSUE_DATA',
          'MODULATOR_COMPLEX', 'TARGET_COMPLEX',
          'MODIFICATIONA', 'MODASEQ', 'MODIFICATIONB', 'MODBSEQ',
          'PMID', 'DIRECT',	'NOTES', 'ANNOTATOR', 'SENTENCE']


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
    if 'dow-reg' in comment:
        breakpoint()
    parts = comment.split(';')
    comment_data = defaultdict(list)
    for part in parts:
        try:
            key, value = part.split(':', maxsplit=1)
        except ValueError:
            print(comment)
            break
        if key not in allowed_keys:
            print(comment)
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
    has_activation = Activation in stmts_by_type
    has_inhibition = Inhibition in stmts_by_type
    has_dephos = Dephosphorylation in stmts_by_type

    if stmt_package[-1].get('mechanism'):
        breakpoint()

    # Here we need to look at cases where there is a SENTENCE
    # involved

    if not (has_activation or has_inhibition):
        if has_dephos:
            for dephos_stmt, dephos_ev, dephos_cur, dephos_comment \
                    in stmts_by_type[Dephosphorylation]:
                if dephos_comment and dephos_comment.get('effect'):
                    # We just handle one effect for now
                    assert len(dephos_comment['effect']) == 1
                    effect = dephos_comment['effect'][0]
                    if 'down-regulates' in effect:
                        stmts_by_type[Inhibition].append(
                            (Inhibition(Agent('X'), Agent('Y')), None, None,
                             {'effect': [effect]})
                        )
                        has_inhibition = True
                        dephos_comment.pop('effect')
                    elif 'up-regulates' in effect:
                        stmts_by_type[Activation].append(
                            (Activation(Agent('X'), Agent('Y')), None, None,
                             {'effect': [effect]})
                        )
                        has_activation = True
                        dephos_comment.pop('effect')
                    else:
                        print('Unknown effect')
                        print(dephos_stmt, dephos_comment)
    elif not has_dephos:
        for act_stmt, act_ev, act_cur, act_comment \
                in stmts_by_type[Activation] + stmts_by_type[Inhibition]:
            if act_comment and act_comment.get('mechanism') == ['dephosphorylation']:
                stmts_by_type[Dephosphorylation].append(
                    (Dephosphorylation(act_stmt.subj, act_stmt.obj), act_ev, act_cur, {})
                )
                has_dephos = True
                act_comment.pop('mechanism')

    if not has_dephos or not (has_activation or has_inhibition):
        return []

    for dephos_stmt_package, activity_stmt_package in \
            itertools.product(stmts_by_type[Dephosphorylation],
                              stmts_by_type[Activation] + stmts_by_type[Inhibition]):
        dephos_stmt, dephos_ev, dephos_cur, dephos_comment = dephos_stmt_package
        activity_stmt, activity_ev, _, activity_comment = activity_stmt_package
        assert not dephos_comment.get('effect')
        phosphatase = dephos_stmt.enz
        substrate = dephos_stmt.sub
        is_activation = isinstance(activity_stmt, Activation)

        if activity_comment and 'effect' in activity_comment:
            assert len(activity_comment['effect']) == 1
            effect = activity_comment['effect'][0]
        else:
            effect = 'up-regulates' if is_activation else 'down-regulates'

        if dephos_stmt.residue and dephos_stmt.position:
            residue = \
                amino_acids[dephos_stmt.residue]['short_name'].capitalize() + \
                dephos_stmt.position
        else:
            residue = ''

        sentence_parts = []
        if dephos_ev.text:
            sentence_parts.append(sanitize_text(dephos_ev.text))
        if activity_ev and activity_ev.text:
            sentence_parts.append(sanitize_text(activity_ev.text))
        sentence = '|'.join(sorted(set(sentence_parts)))

        direct_comment = dephos_comment.get('direct')
        if direct_comment and direct_comment[0].lower() == 'no':
            direct = 'NO'
        else:
            direct = 'YES'

        yield [
            # 'ENTITYA', 'TYPEA', 'IDA', 'DATABASEA'
            phosphatase.name, 'protein', phosphatase.db_refs.get('UP'), 'UNIPROT',
            # 'ENTITYB', 'TYPEB', 'IDB', 'DATABASEB',
            substrate.name, 'protein', substrate.db_refs.get('UP'), 'UNIPROT',
            # 'EFFECT', 'MECHANISM', 'RESIDUE', 'SEQUENCE',
            effect, 'dephosphorylation', residue, '',
            # 'TAX_ID', 'CELL_DATA', 'TISSUE_DATA',
            '9606', '', '',
            # 'MODULATOR_COMPLEX', 'TARGET_COMPLEX',
            '', '',
            # 'MODIFICATIONA', 'MODASEQ', 'MODIFICATIONB', 'MODBSEQ',
            '', '', '', '',
            # 'PMID', 'DIRECT', 'NOTES', 'ANNOTATOR'
            ev.pmid, direct, '', 'lperfetto',
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
    subj = stmt.subj if isinstance(stmt, RegulateActivity) else stmt.enz
    obj = stmt.obj if isinstance(stmt, RegulateActivity) else stmt.sub
    return subj.name, obj.name, ev.pmid


if __name__ == '__main__':
    if os.path.exists('curations.json'):
        with open('curations.json', 'r') as fh:
            curs = json.load(fh)
    else:
        curs = indra_db_rest.get_curations()
        with open('curations.json', 'w') as fh:
            json.dump(curs, fh, indent=1)
    curs = [cur for cur in curs
            if cur.get('source') == 'signor_dephos']
    # Sometimes we have duplicate curations that we can
    # squash here
    curs = {(cur['pa_hash'], cur['source_hash']): cur for cur in curs}.values()
    print('Found %d curations' % len(curs))
    with open('dephosphorylations_with_reg_sorted.pkl', 'rb') as fh:
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

    with open('dephosphorylations_with_reg_export.csv', 'wt') as fh:
        writer = csv.writer(fh, delimiter=',', quotechar='"')
        writer.writerow(header)
        for row in all_rows:
            writer.writerow(row)