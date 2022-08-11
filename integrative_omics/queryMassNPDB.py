#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import os
import sys
import argparse
import sqlite3
import pandas as pd
from multiprocessing import Pool

import gizmos


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Two-column csv with header. One mass signature per line. Column 1: mass signature ID.' 'Column 2: mz.')
    parser.add_argument('output_folder')
    parser.add_argument('--NPDB', default='', required=False, type=str, help='NPDB sqlite file.')
    parser.add_argument('--molecules_table', default='', required=False, type=str, help='CSV with molecules to query. Must have header. Format: name,mm,smiles')
    parser.add_argument('--threads', '-t', default=4, type=int, required=False)
    parser.add_argument('--tolerance', default=30, required=False, type=float, help='In ppm. Default: 30.')
    parser.add_argument('--ion_mode', default='', type=str, required=False, choices=['positive', 'negative'], help='Default: both')
    parser.add_argument('--verbose', '-v', default=False, action='store_true', required=False)
    parser.add_argument('--dev', default=False, action='store_true', required=False, help='Developer mode.')
    return parser.parse_args()

# NPDB = structure_id,isomeric_smiles,canonical_smiles,monoisotopic_mass,molecular_formula,inchi,inchi_key_rdkit,
# NPDB2= Structure_source_id,Structure_id,Source_id,Source_name
# NPDB3= structure_id, monoisotopic_mass, inchi, smile, inchi_key1, inchi_key1, molecular_formula, kingdom, superclass, class, subclass

def NPDB_to_pd(npdb):
    conn = sqlite3.connect(npdb)
    c = conn.cursor()
    query = c.execute("SELECT structure.structure_id,structure.monoisotopic_mass, structure_has_data_source.source_id, structure_has_data_source.source_name, structure.inchi,structure.inchi_key2,structure.smile FROM structure left join structure_has_data_source on structure_has_data_source.structure_id = structure.structure_id")
    cols = [column[0] for column in query.description]
    results = pd.DataFrame.from_records(data=query.fetchall(), columns=cols)
    print(results)
    c.close()
    return results

def query_NPDB(cur_mass_row, npdb):
    # QUERY
    # select * from structure where structure.monoisotopic_mass > 131.96 and structure.monoisotopic_mass < 131.97
    # structure_id,isomeric_smiles,canonical_smiles,monoisotopic_mass,molecular_formula,inchi,inchi_key_rdkit,
    # inchi_key_molconvert,cf_direct_parent,cf_kingdom,cf_superclass,cf_class,cf_subclass,
    # cf_intermediate_0,cf_intermediate_1,cf_intermediate_2,cf_intermediate_3,cf_intermediate_4,cf_intermediate_5,
    # cf_molecular_framework,cf_alternative_parents,cf_substituents,cf_description,cf_queryID
    # mass_range = (cur_mass_row.mm_low, cur_mass_row.mm_high)
    results = []

    for row in npdb.values:
        if (row.monoisotopic_mass >= cur_mass_row.mm_low) and (row.monoisotopic_mass <= cur_mass_row.mm_high):
            temp = [cur_mass_row.ms_name, cur_mass_row.mz, cur_mass_row.observed_mm, cur_mass_row.adduct_name, row[0], row[1], row[2], row[3], row[4], row[5], row[6]]
            print(temp)
            results.append(temp)
        else:
            pass
    labels = ['ms_name', 'mz', 'observed_mm', 'adduct_name', 'NPDB_ID', 'NPDB_mm', 'source_id', 'source_name', 'InChI', 'InChI_key', 'SMILES']
    results = pd.DataFrame(results, columns=labels)
    print(results.columns)
    return results
    # for row in c.execute("SELECT structure.structure_id,structure.monoisotopic_mass, structure_has_data_source.source_id, structure_has_data_source.source_name, structure.inchi,structure.inchi_key2,structure.smile FROM structure left join structure_has_data_source on structure_has_data_source.structure_id = structure.structure_id WHERE structure.monoisotopic_mass >= ? and structure.monoisotopic_mass <= ?", mass_range):
    #    temp = [cur_mass_row.ms_name, cur_mass_row.mz, cur_mass_row.observed_mm, cur_mass_row.adduct_name, row[0], row[1], row[2], row[3], row[4], row[5], row[6]]
    #    print(temp)
    #    results.append(temp)

    # labels = ['ms_name', 'mz', 'observed_mm', 'adduct_name', 'NPDB_ID', 'NPDB_mm', 'source_id', 'source_name', 'InChI', 'InChI_key', 'SMILES']
    # results = pd.DataFrame(results, columns=labels)
    # return results


def query_csv(cur_mass_row, molecules_df):
    mask_low = molecules_df.molecule_mm >= cur_mass_row.mm_low
    mask_high = molecules_df.molecule_mm <= cur_mass_row.mm_high

    molecules_match = molecules_df[mask_low & mask_high].copy()

    if not molecules_match.empty:
        molecules_match['ms_name'] = cur_mass_row.ms_name
        results_df = pd.merge(cur_mass_row.to_frame().transpose(), molecules_match, how='outer')
        results_df = results_df.drop(columns=['mm_low', 'mm_high'])
        return results_df
    else:
        return pd.DataFrame()


#################
# MAIN PIPELINE #
#################
def main():
    if not Options.NPDB and not Options.molecules_table:
        sys.exit('Input the location of the NPDB or a molecules table to query.')

    # OUTPUT FOLDER CREATION & LOG
    gizmos.log_init(Options)

    gizmos.print_milestone('Loading input...', Options.verbose)

    # LOAD IONS FILE
    adducts_df = pd.read_csv(Options.adducts_file, index_col=0)

    if Options.ion_mode:
        adducts_df = adducts_df[adducts_df.Ion_mode == Options.ion_mode]

    # LOAD INPUT
    df = pd.read_csv(Options.input_file, index_col=None, header=0)
    df.columns = ['ms_name', 'mz']
    df = df.astype({'mz': float})

    # GETTING ADDUCTS
    gizmos.print_milestone('Identifying adducts...', Options.verbose)
    if Options.dev:
        df = df.apply(gizmos.get_adduct_data, adducts_df=adducts_df, Options=Options, axis=1)
        df = pd.concat(df.values.tolist()).reset_index(drop=True)
    else:
        results = []
        with Pool(processes=Options.threads) as pool:
            for i, cur_mass_row in df.iterrows():
                pool.apply_async(gizmos.get_adduct_data, args=(cur_mass_row, adducts_df, Options), callback=results.append)
            pool.close()
            pool.join()
        df = pd.concat(results).reset_index(drop=True)
        del results

    # LOAD SQL
    if Options.molecules_table:
        main_cols = ['ms_name', 'molecule_id', 'observed_mm', 'SMILES', 'adduct_name']

        gizmos.print_milestone('Loading input db...', Options.verbose)
        molecules_df = pd.read_csv(Options.molecules_table, header=0, index_col=None)
        molecules_df.rename(columns={molecules_df.columns[0]: 'molecule_id',
                                     molecules_df.columns[1]: 'molecule_mm',
                                     molecules_df.columns[2]: 'SMILES'}, inplace=True)

        gizmos.print_milestone('Querying...', Options.verbose)
        if Options.dev:
            results_df = df.apply(query_csv, molecules_df=molecules_df, axis=1)
            results_df = pd.concat(results_df.values.tolist()).reset_index(drop=True)
        else:
            results_df = []
            with Pool(processes=Options.threads) as pool:
                for i, cur_mass_row in df.iterrows():
                    pool.apply_async(query_csv, args=(cur_mass_row, molecules_df),
                                     callback=results_df.append)
                pool.close()
                pool.join()
            results_df = pd.concat(results_df).reset_index(drop=True)

    else:
        # Options.NPDB
        # Kumar 18th Feb 2022
        # query the NPDB using a specific sql statement and read it in to a pandas dataframe using NPDB_to_pd() function
        # then iterate over the rows of this df for every row of cur_mass_row
        main_cols = ['ms_name', 'NPDB_ID', 'observed_mm', 'SMILES', 'adduct_name', 'InChI']
        gizmos.print_milestone('Connecting to DB and fetching a dataframe..', Options.verbose)
        npdb = NPDB_to_pd(Options.NPDB)
        gizmos.print_milestone('Querying the NPDB DF...', Options.verbose)
        results_df = []
        with Pool(processes=Options.threads) as pool:
            for i, cur_mass_row in df.values:
                pool.apply_async(query_NPDB, args=(cur_mass_row, npdb), callback=results_df.append)
            pool.close()
            pool.join()

        results_df = pd.concat(results_df).reset_index(drop=True)

        if not results_df.empty:
            gizmos.print_milestone('Outputting...', Options.verbose)
            res_file = os.path.join(Options.output_folder, 'full.results.csv')
            results_df.to_csv(res_file, sep=',', index=False)
            # ^^ ms_name, mz, observed_mm, NPDB_ID, NPDB_mm, source_id, source_name, InChI, InChI_key, SMILES

            res_file = os.path.join(Options.output_folder, 'NPDB_entries.txt')
            results_df[['NPDB_ID', 'NPDB_mm', 'SMILES']].drop_duplicates('NPDB_ID').to_csv(res_file, index=False)
            # ^^ NPDB_ID, NPDB_mm, SMILES

    if not results_df.empty:
        # for main results, we only want one entry per ms_name_npdb_id
        res_file = os.path.join(Options.output_folder, 'predicted_structures.csv')
        results_df[main_cols].to_csv(res_file, sep=',', index=False)
        # ^^ ms_name, NPDB_ID, NPDB_mm, SMILES, InChI

        res_file = os.path.join(Options.output_folder, 'observed_mm.csv')
        results_df[['ms_name', 'mz', 'observed_mm']].drop_duplicates().to_csv(res_file, sep=',', index=False)
        # ^^ ms_name, mz, mm

    return


if __name__ == "__main__":
    Options = get_args()
    Options.adducts_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'docs', 'ESI-MS-adducts.csv')
    Options.log_file = os.path.join(Options.output_folder, 'log.txt')
    main()

# python queryMassNPDB.py --NPDB NPDB.sqlite --threads 4 --verbose --dev metabolome.data.csv massNPDB