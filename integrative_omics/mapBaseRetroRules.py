#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import os
import argparse
import pandas as pd
from multiprocessing import Pool
from rdkit import Chem
from rdkit import rdBase

import gizmos

rdBase.DisableLog('rdApp.*')


def get_args():
    """
    Self-explanatory. Toodles.
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',
                        help='File with validated rules.')
    parser.add_argument('output_folder')
    parser.add_argument('--threads', '-t',
                        default=4,
                        type=int,
                        required=False)
    parser.add_argument('--verbose', '-v',
                        default=False,
                        action='store_true',
                        required=False)
    return parser.parse_args()


def query_bigger_smarts(cur_smarts, bigger_df):
    """
    For each cur_smarts it will identify which smarts are bigger, and will query them to see if cur_smarts is in each.
    Results are returned on a ;-separated string, or as "none" (NA).
    :param cur_smarts:
    :param bigger_df:
    :return:
    """
    this_smarts_is_in = bigger_df.mol.apply(lambda x: x.HasSubstructMatch(cur_smarts.mol))  # This is a bool Series
    this_smarts_is_in = bigger_df.smarts_id[this_smarts_is_in]

    if len(this_smarts_is_in):
        this_smarts_is_in = ';'.join([str(n) for n in this_smarts_is_in])
    else:
        this_smarts_is_in = ''
    results = ','.join([cur_smarts.smarts_id, this_smarts_is_in])
    return results


def get_representative_smarts(cur_rule, has_identical_smarts):
    """
    Groups cur_smarts and its identities to elect a representative smarts (smallest ID)
    :param cur_rule:
    :param has_identical_smarts:
    :return:
    """
    if cur_rule.smarts_id in has_identical_smarts:
        smarts_group = list(has_identical_smarts[cur_rule.smarts_id])
        smarts_group.append(cur_rule.smarts_id)
        smarts_group = sorted(smarts_group)
        return smarts_group[0]
    else:
        return cur_rule.smarts_id


def find_base_rules(cur_rule, has_smaller_smarts, has_identical_smarts):
    # ## Case 1 = yes has_smaller &             yes identical   =v
    # ## ## Case 1.1 = smaller == identity => choose representative base rule (smallest smarts_id)
    # ## ## Case 1.2 = smaller ~= identity => NOT base rule
    # ## Case 2 = yes has_smaller &             no  identical   => NOT base rule
    # ## Case 3 = no  has_smaller & yes is_in & yes identical   => choose representative base rule (smallest smarts_id)
    # ## Case 4 = no  has_smaller & yes is_in & no  identical   => base rule
    # ## Case 5 = no  has_smaller & no  is_in & yes identical   => choose representative base rule (smallest smarts_id)
    # ## Case 6 = no  has_smaller & no  is_in & no  identical   => base rule
    cur_rule['representative_smarts'] = get_representative_smarts(cur_rule, has_identical_smarts)

    if cur_rule.smarts_id in has_smaller_smarts:
        cur_rule['smarts_has'] = has_smaller_smarts[cur_rule.smarts_id]
        if cur_rule.smarts_id in has_identical_smarts:          # Case 1
            cur_rule['identity'] = has_identical_smarts[cur_rule.smarts_id]
            if has_smaller_smarts[cur_rule.smarts_id] == \
                    has_identical_smarts[cur_rule.smarts_id]:  # Case 1.1
                if cur_rule.smarts_id == cur_rule['representative_smarts']:
                    cur_rule['is_base'] = True
                else:
                    cur_rule['is_base'] = False
            else:                                                                                   # Case 1.2
                cur_rule['is_base'] = False
        else:                                                   # Case 2
            cur_rule['identity'] = None
            cur_rule['is_base'] = False
    elif cur_rule.smarts_id in has_identical_smarts:            # Case 3 & 5
        cur_rule['smarts_has'] = None
        cur_rule['identity'] = has_identical_smarts[cur_rule.smarts_id]
        if cur_rule.smarts_id == cur_rule['representative_smarts']:
            cur_rule['is_base'] = True
        else:
            cur_rule['is_base'] = False
    else:                                                       # Case 4 & 6
        cur_rule['smarts_has'] = None
        cur_rule['identity'] = None
        cur_rule['is_base'] = True

    # clean sets so only identity has identity and then convert to strings
    cur_rule['smarts_is_in'] = gizmos.pd_to_set(cur_rule.smarts_is_in, ';')
    if not pd.isnull(cur_rule.identity):
        if not pd.isnull(cur_rule.smarts_is_in):
            cur_rule['smarts_is_in'] = cur_rule.smarts_is_in - cur_rule.identity
        if not pd.isnull(cur_rule.smarts_has):
            cur_rule['smarts_has'] = cur_rule.smarts_has - cur_rule.identity

    # convert to strings
    cur_rule['identity'] = gizmos.set_to_string(cur_rule.identity, ';')
    cur_rule['smarts_has'] = gizmos.set_to_string(cur_rule.smarts_has, ';')
    cur_rule['smarts_is_in'] = gizmos.set_to_string(cur_rule.smarts_is_in, ';')

    return cur_rule


def write_output(results_row):
    """
    Writes output.
    :param results_row:
    :return:
    """
    res_file = os.path.join(Options.output_folder, 'all_small_rules.csv')
    with open(res_file, 'a') as resfile:
        resfile.write(results_row)
        resfile.write('\n')

    if not results_row.endswith(','):       # If there is a smarts_is_in
        network_file = os.path.join(Options.output_folder, 'network.csv')
        sub_smarts = results_row.split(',')[0]
        with open(network_file, 'a') as networkfile:
            for super_smarts in results_row.split(',')[-1].split(';'):
                line = ','.join([sub_smarts, super_smarts])
                networkfile.write(line)
                networkfile.write('\n')
    return


#################
# MAIN PIPELINE #
#################
def main():
    # OUTPUT FOLDER CREATION & LOG
    gizmos.log_init(Options)

    # LOAD INPUT
    gizmos.print_milestone('Loading input...', Options.verbose)
    df = pd.read_csv(Options.input_file, index_col=None,
                     dtype={'reaction_id': str, 'substrate_id': str, 'validated': bool})
    # ^^ reaction_id,substrate_id,diameter,direction,smarts_id,rxn_smarts,validated
    df = df[df.validated].reset_index(drop=True)                        # Ignore unvalidated rules
    df['reaction_substrate'] = df.reaction_id + '_' + df.substrate_id

    # DEFINING SMALL RULES
    gizmos.print_milestone('Defining small rules...', Options.verbose)
    all_small_rules = []
    for cur_reaction_substrate in df.reaction_substrate.unique():       # Loop through reaction_substrate
        mask = df.reaction_substrate == cur_reaction_substrate
        index_of_smallest_diameter = df[mask].sort_values(by='diameter').iloc[0].name   # sort ascending by diameter
        all_small_rules.append(index_of_smallest_diameter)                              # ..and pick first row.

    all_small_rules = df.loc[all_small_rules].copy().reset_index(drop=True)        # This is a df of small_rules only
    del df

    gizmos.print_milestone(str(len(all_small_rules)) + ' structures to query found...', Options.verbose)

    # PRE-PROCESSING
    gizmos.print_milestone("Mol'ing...", Options.verbose)
    all_small_rules['substrate_str'] = all_small_rules.rxn_smarts.apply(gizmos.get_subs_string)
    # Create a df specifically for mols so we don't mol same structure twice. Index=smarts_id
    mols_df = all_small_rules.drop_duplicates(subset=['smarts_id'], keep='first')[['smarts_id', 'substrate_str']]
    mols_df = mols_df.set_index('smarts_id')
    mols_df['mol'] = mols_df.substrate_str.apply(Chem.MolFromSmarts)
    mols_df['mm'] = mols_df.mol.apply(gizmos.get_mm_from_mol, is_smarts=True)   # Mass of the pattern.
    # Merge mols_df into original df to continue normal pipeline.
    del mols_df['substrate_str']    # to avoid duplicate column
    all_small_rules = pd.merge(all_small_rules, mols_df, left_on='smarts_id', right_index=True)
    del mols_df

    # INITIALIZING OUTPUT
    res_file = os.path.join(Options.output_folder, 'all_small_rules.csv')
    with open(res_file, 'w') as resfile:
        resfile.write('smarts_id,smarts_is_in')
        resfile.write('\n')
    network_file = os.path.join(Options.output_folder, 'network.csv')
    with open(network_file, 'w') as networkfile:
        networkfile.write('subs,super')
        networkfile.write('\n')

    # FINDING SMALL RULES RELATIONSHIPS (IS_IN)
    gizmos.print_milestone('Querying...', Options.verbose)
    with Pool(processes=Options.threads) as pool:
        for i, cur_rule in all_small_rules.iterrows():
            mm_mask = all_small_rules.mm >= cur_rule.mm
            not_same_mask = all_small_rules.smarts_id != cur_rule.smarts_id
            bigger_df = all_small_rules[mm_mask & not_same_mask]
            pool.apply_async(query_bigger_smarts, args=(cur_rule, bigger_df), callback=write_output)
        pool.close()
        pool.join()
    del all_small_rules

    # DEFINING BASE RULES
    small_rules_df = pd.read_csv(res_file, index_col=None, dtype={'smarts_id': str, 'smarts_is_in': str})
    small_rules_df = small_rules_df.drop_duplicates(subset='smarts_id').set_index('smarts_id')
    small_rules_df.to_csv(res_file, index=True)

    # # find identities
    gizmos.print_milestone('Processing identities...', Options.verbose)
    has_smaller_smarts = {}
    has_identical_smarts = {}
    for smarts_id, cur_rule in small_rules_df.iterrows():
        if pd.isnull(cur_rule.smarts_is_in):                                    # no matches = no identities
            continue
        else:
            for bigger_smarts in cur_rule.smarts_is_in.split(';'):                  # loop through all bigger_smarts
                if bigger_smarts not in has_smaller_smarts:
                    has_smaller_smarts[bigger_smarts] = set()                       # Init
                has_smaller_smarts[bigger_smarts].add(smarts_id)                     # Here we fill smarts_has
                if (not pd.isnull(small_rules_df.smarts_is_in[bigger_smarts]) and    # bigger_smarts has bigger_smarts?
                        smarts_id in small_rules_df.smarts_is_in[bigger_smarts]):  # bigger_smarts is cur_smarts?
                    if bigger_smarts not in has_identical_smarts:
                        has_identical_smarts[bigger_smarts] = set()
                    if smarts_id not in has_identical_smarts:
                        has_identical_smarts[smarts_id] = set()
                    has_identical_smarts[bigger_smarts].add(smarts_id)              # then we add to identities
                    has_identical_smarts[smarts_id].add(bigger_smarts)

    # # find base rules
    # ## Case 1 = yes has_smaller &             yes identical   =v
    # ## ## Case 1.1 = smaller == identity => choose representative base rule (smallest smarts_id)
    # ## ## Case 1.2 = smaller ~= identity => NOT base rule
    # ## Case 2 = yes has_smaller &             no  identical   => NOT base rule
    # ## Case 3 = no  has_smaller & yes is_in & yes identical   => choose representative base rule (smallest smarts_id)
    # ## Case 4 = no  has_smaller & yes is_in & no  identical   => base rule
    # ## Case 5 = no  has_smaller & no  is_in & yes identical   => choose representative base rule (smallest smarts_id)
    # ## Case 6 = no  has_smaller & no  is_in & no  identical   => base rule
    gizmos.print_milestone('Finding base rules...', Options.verbose)
    small_rules_df = small_rules_df.rename_axis('smarts_id').reset_index()
    small_rules_df = small_rules_df.apply(find_base_rules,
                                          args=(has_smaller_smarts, has_identical_smarts,), axis=1)

    small_rules_df = small_rules_df[['smarts_id', 'smarts_is_in', 'smarts_has', 'identity',
                                     'representative_smarts', 'is_base']]
    res_file = os.path.join(Options.output_folder, 'base_rules.csv')
    small_rules_df.to_csv(res_file, index=False)
    # ^^ smarts_id, smarts_is_in, smarts_has, identity, representative_smarts, is_base

    return


if __name__ == "__main__":
    Options = get_args()
    Options.log_file = os.path.join(Options.output_folder, 'log.txt')
    main()
