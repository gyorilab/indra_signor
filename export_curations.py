import csv
import pickle
import itertools
from collections import defaultdict
from indra.statements import Dephosphorylation, RegulateActivity, Activation, Inhibition
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
    text = text.replace('[XREF_BIBR - XREF_BIBR]', '')
    text = text.replace('[XREF_BIBR]', '')
    text = text.replace('[XREF_BIBR', '')
    text = text.replace('XREF_BIBR', '')
    text = text.strip()
    return text


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

    if not has_dephos or not (has_activation or has_inhibition):
        return []

    for dephos_stmt_package, activity_stmt_package in \
            itertools.product(stmts_by_type[Dephosphorylation],
                              stmts_by_type[Activation] + stmts_by_type[Inhibition]):
        dephos_stmt, dephos_ev, dephos_cur = dephos_stmt_package
        activity_stmt, activity_ev, activity_cur = activity_stmt_package
        phosphatase = dephos_stmt.enz
        substrate = dephos_stmt.sub
        is_activation = isinstance(activity_stmt, Activation)

        effect = 'up-regulates' if is_activation else 'down-regulates'

        if dephos_stmt.residue and dephos_stmt.position:
            residue = \
                amino_acids[dephos_stmt.residue]['short_name'].capitalize() + \
                dephos_stmt.position
        else:
            residue = ''

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
            ev.pmid, 'YES', '', 'lperfetto',
            # 'SENTENCE'
            '|'.join([sanitize_text(dephos_ev.text),
                      sanitize_text(activity_ev.text)]),
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
                    row_wit_site_sentences = set(row_with_site['SENTENCE'].split('|'))
                    row_sentences = set(row['SENTENCE'].split('|'))
                    unique_sentences = '|'.join(list(row_sentences |
                                                     row_wit_site_sentences))
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
    curs = indra_db_rest.get_curations()
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
        stmt = stmts_by_hash[cur['pa_hash']]
        ev = ev_by_source_hash[cur['source_hash']]
        curs_by_hashes[(stmt.get_hash(), ev.get_source_hash())].append(cur)
        key = get_pair_key(stmt, ev)
        stmts_by_pair_pubmed_key[key].append((stmt, ev, cur))

    all_rows = []
    for key, curations in stmts_by_pair_pubmed_key.items():
        print(key)
        for s, ev, cur in curations:
            print(s, ev.text)
        curation_rows = list(curations_to_rows(curations))
        curation_rows = merge_curation_rows(curation_rows)
        all_rows += curation_rows

    with open('dephosphorylations_with_reg_export.csv', 'wt') as fh:
        writer = csv.writer(fh, delimiter=',', quotechar='"')
        writer.writerow(header)
        for row in all_rows:
            writer.writerow(row)