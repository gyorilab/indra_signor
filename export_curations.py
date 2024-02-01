import csv
import pickle
import itertools
from collections import defaultdict
from indra.statements import Dephosphorylation, RegulateActivity, Activation, Inhibition
from indra.statements import amino_acids
from indra.sources import indra_db_rest


header = ['ENTITYA', 'TYPEA', 'IDA', 'DATABASEA'
          'ENTITYB', 'TYPEB', 'IDB', 'DATABASEB',
          'EFFECT', 'MECHANISM', 'RESIDUE', 'SEQUENCE',
          'TAX_ID', 'CELL_DATA', 'TISSUE_DATA',
          'MODULATOR_COMPLEX', 'TARGET_COMPLEX',
          'MODIFICATIONA', 'MODASEQ', 'MODIFICATIONB', 'MODBSEQ',
          'PMID', 'DIRECT',	'NOTES', 'ANNOTATOR', 'SENTENCE']


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
            # 'PMID', 'DIRECT', 'NOTES', 'ANNOTATOR', 'SENTENCE',
            ev.pmid, 'YES', '', 'lperfetto', dephos_ev.text + activity_ev.text,
        ]


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
        all_rows += curations_to_rows(curations)

    with open('dephosphorylations_with_reg_export.csv', 'wt') as fh:
        writer = csv.writer(fh, delimiter=',', quotechar='"')
        writer.writerow(header)
        for row in all_rows:
            writer.writerow(row)