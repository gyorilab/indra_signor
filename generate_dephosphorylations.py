import gzip
import csv
import json
import tqdm
import codecs
import pystow
from collections import defaultdict
from indra.statements import stmt_from_json, Dephosphorylation, RegulateActivity
from indra.databases.hgnc_client import is_phosphatase
from indra.util import batch_iter
from indra.sources import indra_db_rest
import indra.tools.assemble_corpus as ac
from indra.assemblers.html import HtmlAssembler


processed_stmts_path = pystow.join('indra', 'db',
                                   name='processed_statements.tsv.gz')

def load_statement_json(json_str: str, attempt: int = 1, max_attempts: int = 5):
    try:
        return json.loads(json_str)
    except json.JSONDecodeError:
        if attempt < max_attempts:
            json_str = codecs.escape_decode(json_str)[0].decode()
            return load_statement_json(
                json_str, attempt=attempt + 1, max_attempts=max_attempts
            )
    raise ValueError

# This function can be used to get a clean set of just
# dephosphorylation statements but since SIGNOR
# requires polarity of regulation along with PTMs we need
# to use get_dephos_act_inh_stmts instead
def get_dephos_stmts():
    # Initial filter for dephosphorylations
    stmts_by_hash = defaultdict(list)
    with gzip.open(processed_stmts_path, 'rt') as fh:
        reader = csv.reader(fh, delimiter='\t')
        for lines in tqdm.tqdm(batch_iter(reader, 10000), total=6041):
            for stmt_hash, stmt_json in lines:
                if 'dephosphorylation' not in stmt_json.lower():
                    continue
                stmt_json = load_statement_json(stmt_json)
                stmt = stmt_from_json(stmt_json)
                if isinstance(stmt, Dephosphorylation):
                    if stmt.enz is not None:
                        if is_phosphatase(stmt.enz.name):
                            stmts_by_hash[stmt_hash].append(stmt)
    return stmts_by_hash


def get_dephos_act_inh_stmts():
    stmts_by_hash = defaultdict(list)
    with gzip.open(processed_stmts_path, 'rt') as fh:
        reader = csv.reader(fh, delimiter='\t')
        for lines in tqdm.tqdm(batch_iter(reader, 10000), total=6041):
            for stmt_hash, stmt_json in lines:
                stmt_json = load_statement_json(stmt_json)
                stmt = stmt_from_json(stmt_json)
                if isinstance(stmt, (Dephosphorylation, RegulateActivity)):
                    agents = stmt.agent_list()
                    if len(agents) != 2 or agents[0] is None:
                        continue
                    controller = agents[0]
                    if controller.db_refs.get('HGNC'):
                        if is_phosphatase(controller.name):
                            stmts_by_hash[stmt_hash].append(stmt)
    return stmts_by_hash


if __name__ == '__main__':
    stmts_by_hash = get_dephos_act_inh_stmts()
    # Filter for human genes, non-self regulations
    stmts_by_hash = {k: v for k, v in stmts_by_hash.items() if
                     'HGNC' in v[0].real_agent_list()[1].db_refs}
    stmts_by_hash = {k: v for k, v in stmts_by_hash.items() if
                     v[0].real_agent_list()[0].name != v[0].real_agent_list()[1].name}
    # Flatten list of hash-based lists
    stmts = []
    for ss in stmts_by_hash.values():
        stmts += ss

    # Figure out the sets of Statements and evidences where
    # a dephosphorylation statement and an activation/inhibition
    # statement come from the same paper
    stmts_by_paper = defaultdict(list)
    for stmt in stmts:
        for ev in stmt.evidence:
            if ev.pmid:
                if stmt not in stmts_by_paper[ev.pmid]:
                    stmts_by_paper[ev.pmid].append(stmt)

    relevant_stmts_pmids = defaultdict(set)
    for paper, paper_stmts in stmts_by_paper.items():
        stmt_types_by_pair = defaultdict(set)
        for stmt in paper_stmts:
            agents = stmt.real_agent_list()
            stmt_types_by_pair[(agents[0].name, agents[1].name)].add(type(stmt))
        for pair, stmt_types in stmt_types_by_pair.items():
            if Dephosphorylation in stmt_types and len(stmt_types) > 1:
                relevant_stmts_pmids[pair].add(paper)

    # Get just dephosphorylation statements
    dephos_stmts = [stmt for stmt in stmts if isinstance(stmt, Dephosphorylation)]

    # Organize dephos Statements by enzyme/substrate pair
    # but only keep those Statements that contain an Activation/Inhibition
    # from the same paper as Dephos
    dephos_stmts_by_pair = defaultdict(list)
    for stmt in dephos_stmts:
        agents = stmt.real_agent_list()
        pair = (agents[0].name, agents[1].name)
        if pair in relevant_stmts_pmids:
            if stmt.evidence[0].pmid in relevant_stmts_pmids[pair]:
                dephos_stmts_by_pair[pair].append(stmt)

    # Organize SIGNOR Dephos Statements by whether they have a specific
    # site or not
    has_signors_without_site = {}
    has_signors_with_site = defaultdict(set)
    for k, v in dephos_stmts_by_pair.items():
        for stmt in v:
            if stmt.evidence[0].source_api == 'signor':
                if stmt.position:
                    has_signors_with_site[k].add(stmt.position)
                else:
                    has_signors_without_site[k] = True

    # Filter out Statements that are already in SIGNOR
    dephos_stmts_by_pair_new = defaultdict(list)
    for k, v in dephos_stmts_by_pair.items():
        for stmt in v:
            if k in has_signors_with_site:
                if stmt.position is None:
                    continue
                elif stmt.position in has_signors_with_site[k]:
                    continue
            elif k in has_signors_without_site:
                if stmt.position is None:
                    continue
            dephos_stmts_by_pair_new[k].append(stmt)

    # We now compile all relevant Statements based on new dephosphorylations
    all_relevant_stmts = []
    for stmt in stmts:
        agents = stmt.real_agent_list()
        pair = (agents[0].name, agents[1].name)
        if pair in dephos_stmts_by_pair_new:
            if stmt.evidence[0].pmid and \
                    stmt.evidence[0].pmid in relevant_stmts_pmids[pair]:
                if not (isinstance(stmt, Dephosphorylation) and
                        stmt.evidence[0].source_api == 'signor'):
                    all_relevant_stmts.append(stmt)

    # Assemble Statements and dump the results
    stmts_assembled = ac.run_preassembly(all_relevant_stmts,
                                         return_toplevel=False)

    # Filter for curated statements
    curs = indra_db_rest.get_curations()
    stmts_correct = ac.filter_by_curation(stmts_assembled, curs)

    # Filter out incorrect Statements
    #correct_hashes = {c['pa_hash'] for c in curs if c['tag'] == 'correct'}
    #stmts_correct = [s for s in stmts_assembled if s.get_hash() in correct_hashes]

    for stmt in stmts_correct:
        stmt.evidence = sorted(stmt.evidence, key=lambda x: x.pmid)

    # This function sorts statements such that they are grouped according
    # to phosphatase/target pairs, sorted according to the number of
    # evidences for dephosphorylation overall, and then sorted according
    # to the number of evidences for the specific statement type with
    # Dephosphorylation always coming first within the group.
    def get_sort_key(stmt, ev_counts):
        agents = stmt.real_agent_list()
        a, b = [a.name for a in agents[:2]]
        ab_key = (a, b)
        group_ev_key = ev_counts.get(ab_key, 0)
        stmt_type_key = 0 if isinstance(stmt, Dephosphorylation) else 1
        ev_key = len(stmt.evidence)
        return (-group_ev_key, ab_key, stmt_type_key, -ev_key)

    # We need these to sort the statements, these are evidences
    # for dephosphorylation specifically
    ev_counts = defaultdict(int)
    for stmt in stmts_correct:
        agents = stmt.real_agent_list()
        a, b = [a.name for a in agents[:2]]
        ab_key = (a, b)
        if isinstance(stmt, Dephosphorylation):
            ev_counts[ab_key] += len(stmt.evidence)
    ev_counts = dict(ev_counts)

    sorted_stmts = sorted(stmts_correct, key=lambda x: get_sort_key(x, ev_counts))

    # Dump Statements as pickle as HTML
    ac.dump_statements(sorted_stmts, 'dephosphorylations_with_reg_sorted.pkl')
    ha = HtmlAssembler(sorted_stmts,
                       db_rest_url='https://db.indra.bio',
                       title='INDRA Dephosphorylation Statements with regulation sample')
    m = ha.make_model(grouping_level='agent-pair')
    ha.save_model('indra_dephosphorylation_with_reg_sample.html')
