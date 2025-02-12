import os
import gzip
import csv
import json
import tqdm
import codecs
import pystow
from collections import defaultdict
from indra.statements import stmt_from_json, RegulateActivity, RegulateAmount,\
    modtype_to_modclass, Modification
from indra.databases import hgnc_client
from indra.util import batch_iter
from indra.sources import indra_db_rest
import indra.tools.assemble_corpus as ac
from indra.assemblers.html import HtmlAssembler


processed_stmts_path = pystow.join('indra', 'db',
                                   name='processed_statements.tsv.gz')


def enzyme_is_phosphatase(stmt):
    if isinstance(stmt, Modification) and stmt.enz is not None:
        if hgnc_client.is_phosphatase(stmt.enz.name):
            return True
    elif isinstance(stmt, (RegulateActivity, RegulateAmount)):
        if hgnc_client.is_phosphatase(stmt.subj.name):
            return True
    return False


def enzyme_is_kinase(stmt):
    if isinstance(stmt, Modification) and stmt.enz is not None:
        if hgnc_client.is_kinase(stmt.enz.name):
            return True
    elif isinstance(stmt, (RegulateActivity, RegulateAmount)):
        if hgnc_client.is_kinase(stmt.subj.name):
            return True
    return False


def get_e3_ubi_ligases():
    if not os.path.exists('e3_ubi_ligases.txt'):
        df = pystow.ensure_excel(
            'indra_signor',
            url='https://esbl.nhlbi.nih.gov/Databases/KSBP2/'
                'Targets/Lists/E3-ligases/'
                'Definite%20Ligase%20List.xlsx',
            read_excel_kwargs={'skiprows': 12})
        gene_names = []
        for gene in df['Gene Symbol']:
            if hgnc_client.get_hgnc_id(gene):
                gene_names.append(gene)
            else:
                # Some outdated gene names and ones that aren't all caps
                hgnc_id = hgnc_client.get_current_hgnc_id(gene.upper())
                if hgnc_id:
                    current_name = hgnc_client.get_hgnc_name(hgnc_id)
                    gene_names.append(current_name)
                else:
                    print(f'Gene {gene} not found in HGNC')
        with open('e3_ubi_ligases.txt', 'w') as fh:
            for gene_name in gene_names:
                fh.write(f'{gene_name}\n')
    else:
        with open('e3_ubi_ligases.txt', 'r') as fh:
            gene_names = [line.strip()
                          for line in fh.readlines() if line.strip()]
    return gene_names


UBI_LIGASES = get_e3_ubi_ligases()


def enzyme_is_e3_ubi_ligase(stmt):
    assert len(stmt.evidence) == 1
    if stmt.evidence[0].text == 'UbiNet 2.0':
        return False
    if isinstance(stmt, Modification) and stmt.enz is not None:
        if stmt.enz.name in UBI_LIGASES:
            return True
    elif isinstance(stmt, (RegulateActivity, RegulateAmount)):
        if stmt.subj.name in UBI_LIGASES:
            return True
    return False


semantic_constraints = {
    'dephosphorylation': enzyme_is_phosphatase,
    'phosphorylation': enzyme_is_kinase,
    'ubiquitination': enzyme_is_e3_ubi_ligase
}


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
# a given type of modification statements but since SIGNOR
# requires polarity of regulation along with PTMs we need
# to use get_mod_and_reg_stmts instead
def get_mod_stmts(mod_key='dephophorylation'):
    # Initial filter for a given mod
    stmts_by_hash = defaultdict(list)
    with gzip.open(processed_stmts_path, 'rt') as fh:
        reader = csv.reader(fh, delimiter='\t')
        for lines in tqdm.tqdm(batch_iter(reader, 10000), total=6041):
            for stmt_hash, stmt_json in lines:
                if mod_key not in stmt_json.lower():
                    continue
                stmt_json = load_statement_json(stmt_json)
                stmt = stmt_from_json(stmt_json)
                if isinstance(stmt, modtype_to_modclass[mod_key]):
                    constraint = semantic_constraints.get(mod_key)
                    if constraint:
                        if constraint(stmt):
                            stmts_by_hash[stmt_hash].append(stmt)
                    else:
                        stmts_by_hash[stmt_hash].append(stmt)
    return stmts_by_hash


def get_mod_and_reg_stmts(mod_key='dephophorylation', regulation_types=None):
    if regulation_types is None:
        regulation_types = {RegulateActivity}
    stmts_by_hash = defaultdict(list)
    stmt_types_allowed = \
        tuple({modtype_to_modclass[mod_key]} | regulation_types)
    with gzip.open(processed_stmts_path, 'rt') as fh:
        reader = csv.reader(fh, delimiter='\t')
        for lines in tqdm.tqdm(batch_iter(reader, 10000), total=6041):
            for stmt_hash, stmt_json in lines:
                stmt_json = load_statement_json(stmt_json)
                stmt = stmt_from_json(stmt_json)
                if isinstance(stmt, stmt_types_allowed):
                    agents = stmt.agent_list()
                    if len(agents) != 2 or agents[0] is None:
                        continue
                    controller = agents[0]
                    if controller.db_refs.get('HGNC'):
                        constraint = semantic_constraints.get(mod_key)
                        if constraint:
                            if constraint(stmt):
                                stmts_by_hash[stmt_hash].append(stmt)
                        else:
                            stmts_by_hash[stmt_hash].append(stmt)
    return stmts_by_hash


if __name__ == '__main__':
    ####################
    mod_key = 'phosphorylation'
    regulation_types = {RegulateAmount, RegulateActivity}
    ####################
    mod_class = modtype_to_modclass[mod_key]
    stmts_by_hash = get_mod_and_reg_stmts(mod_key, regulation_types)
    print(f'Got {len(stmts_by_hash)} statements of type {mod_key} and'
          f' regulation types {regulation_types}')
    # Filter for human genes, non-self regulations
    stmts_by_hash = {k: v for k, v in stmts_by_hash.items() if
                     'HGNC' in v[0].real_agent_list()[1].db_refs}
    print(f'Filtered to {len(stmts_by_hash)} statements with human genes')
    stmts_by_hash = {k: v for k, v in stmts_by_hash.items() if
                     v[0].real_agent_list()[0].name != v[0].real_agent_list()[1].name}
    print(f'Filtered to {len(stmts_by_hash)} statements with non-self regulations')
    # Flatten list of hash-based lists
    stmts = []
    for ss in stmts_by_hash.values():
        stmts += ss

    # Filter for curated statements
    curs = indra_db_rest.get_curations()
    stmts = ac.filter_by_curation(stmts, curs)

    # Get just modification statements of the given type
    mod_stmts = [stmt for stmt in stmts if isinstance(stmt, mod_class)]

    # Organize mod Statements by enzyme/substrate pair
    mod_stmts_by_pair = defaultdict(list)
    for stmt in mod_stmts:
        agents = stmt.real_agent_list()
        pair = (agents[0].name, agents[1].name)
        mod_stmts_by_pair[pair].append(stmt)

    # Organize SIGNOR modification Statements by whether they have a
    # specific site or not
    has_signors_without_site = {}
    has_signors_with_site = defaultdict(set)
    for k, v in mod_stmts_by_pair.items():
        for stmt in v:
            if stmt.evidence[0].source_api == 'signor':
                if stmt.position:
                    has_signors_with_site[k].add(stmt.position)
                else:
                    has_signors_without_site[k] = True
    has_signors_with_site = dict(has_signors_with_site)

    # Filter out Statements that are already in SIGNOR
    mod_stmts_by_pair_new = defaultdict(set)
    for k, v in mod_stmts_by_pair.items():
        for stmt in v:
            if k in has_signors_with_site:
                if stmt.position is None:
                    continue
                elif stmt.position in has_signors_with_site[k]:
                    continue
            elif k in has_signors_without_site:
                if stmt.position is None:
                    continue
            mod_stmts_by_pair_new[k].add(stmt)
    mod_stmts_by_pair_new = dict(mod_stmts_by_pair_new)

    # We now organize all statements by relevant paper
    stmts_by_paper = defaultdict(dict)
    for stmt in stmts:
        pair = (stmt.real_agent_list()[0].name, stmt.real_agent_list()[1].name)
        if pair not in mod_stmts_by_pair_new:
            continue
        for ev in stmt.evidence:
            if not ev.pmid:
                continue
            if pair in stmts_by_paper[ev.pmid]:
                stmts_by_paper[ev.pmid][pair].add(stmt)
            else:
                stmts_by_paper[ev.pmid][pair] = {stmt}
    stmts_by_paper = dict(stmts_by_paper)

    all_relevant_stmts = []
    for paper, stmts_by_pair in stmts_by_paper.items():
        for pair, stmts_for_pair in stmts_by_pair.items():
            stmts_by_type = defaultdict(list)
            for stmt in stmts_for_pair:
                if isinstance(stmt, mod_class):
                    if stmt in mod_stmts_by_pair_new[pair]:
                        stmts_by_type[type(stmt)].append(stmt)
                else:
                    stmts_by_type[type(stmt)].append(stmt)
            # This means that in this paper, for this pair, we have
            # both a modification and a regulation statement
            if mod_class in stmts_by_type and len(stmts_by_type) > 1:
                # Now we need to save all _novel_ modifications and all
                # regulation statements from this paper
                for stmt in stmts_for_pair:
                    if isinstance(stmt, mod_class):
                        if stmt in mod_stmts_by_pair_new[pair]:
                            all_relevant_stmts.append(stmt)
                    else:
                        all_relevant_stmts.append(stmt)

    # Assemble Statements and dump the results
    stmts_assembled = ac.run_preassembly(all_relevant_stmts,
                                         return_toplevel=False)

    # Sort evidence by PMID
    for stmt in stmts_assembled:
        stmt.evidence = sorted(stmt.evidence, key=lambda x: x.pmid)

    # This function sorts statements such that they are grouped according
    # to phosphatase/target pairs, sorted according to the number of
    # evidences for the given type of modification overall, and then sorted
    # according to the number of evidences for the specific statement
    # type with Modification always coming first within the group.
    def get_sort_key(stmt, ev_counts):
        agents = stmt.real_agent_list()
        a, b = [a.name for a in agents[:2]]
        ab_key = (a, b)
        group_ev_key = ev_counts.get(ab_key, 0)
        stmt_type_key = 0 if isinstance(stmt, mod_class) else 1
        ev_key = len(stmt.evidence)
        return (-group_ev_key, ab_key, stmt_type_key, -ev_key)

    # We need these to sort the statements, these are evidences
    # for the given type of modification specifically
    ev_counts = defaultdict(int)
    for stmt in stmts_assembled:
        agents = stmt.real_agent_list()
        a, b = [a.name for a in agents[:2]]
        ab_key = (a, b)
        if isinstance(stmt, mod_class):
            ev_counts[ab_key] += len(stmt.evidence)
    ev_counts = dict(ev_counts)

    sorted_stmts = sorted(stmts_assembled, key=lambda x: get_sort_key(x, ev_counts))

    # Dump Statements as pickle as HTML
    ac.dump_statements(sorted_stmts, f'{mod_key}s_with_reg_sorted.pkl')
    ha = HtmlAssembler(sorted_stmts,
                       db_rest_url='https://db.indra.bio',
                       title=f'INDRA {mod_key} Statements with regulation sample')
    m = ha.make_model(grouping_level='agent-pair')
    ha.save_model(f'indra_{mod_key}_with_reg.html')
