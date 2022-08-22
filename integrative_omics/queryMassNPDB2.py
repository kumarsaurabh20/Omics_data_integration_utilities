import io
import os
import time
import sqlite3
import pandas as pd
import numpy as np
from multiprocessing import Pool


def query_NPDB(cur_mass_row, npdb):
    results = []

    for row in npdb.to_dict('records'):
        if (row['monoisotopic_mass'] >= cur_mass_row[3]) and (row['monoisotopic_mass'] <= cur_mass_row[4]):
            temp = [cur_mass_row[0], cur_mass_row[1], cur_mass_row[2], cur_mass_row[5], row['structure_id'], row['monoisotopic_mass'], row['source_id'], row['source_name'], row['inchi'], row['inchi_key2'], row['smile']]
            results.append(temp)
        else:
            pass
    labels = ['ms_name', 'mz', 'observed_mm', 'adduct_name', 'NPDB_ID', 'NPDB_mm', 'source_id', 'source_name', 'InChI', 'InChI_key', 'SMILES']
    results = pd.DataFrame(results, columns=labels)
    #print(results.columns)
    return results


def NPDB_to_pd(npdb):
    conn = sqlite3.connect(npdb)
    c = conn.cursor()
    query = c.execute("SELECT structure.structure_id,structure.monoisotopic_mass, structure_has_data_source.source_id, structure_has_data_source.source_name, structure.inchi,structure.inchi_key2,structure.smile FROM structure left join structure_has_data_source on structure_has_data_source.structure_id = structure.structure_id")
    cols = [column[0] for column in query.description]
    results = pd.DataFrame.from_records(data=query.fetchall(), columns=cols)
    c.close()
    temp = results.head(n=100)
    temp.to_csv("table4.csv")
    return results


def get_mass_range(cur_mass_row, tolerance):

    mz = cur_mass_row[1]
    tol = mz * (tolerance / 1000000)
    top = mz + tol
    btm = mz - tol
    temp = np.append(cur_mass_row, [btm, top], axis=0)
    return temp


def get_adduct_data(cur_mass_row, adducts_df, tolerance):
    """
    Gets possible mm's for the mz of cur_mass
    """
    def get_cur_mm_adducts(cur_adduct_row):
        mm_low = (cur_mass_row[2] - cur_adduct_row[5]) / cur_adduct_row[4]
        mm_high = (cur_mass_row[3] - cur_adduct_row[5]) / cur_adduct_row[4]
        mm = (cur_mass_row[1] - cur_adduct_row[5]) / cur_adduct_row[4]
        adduct_name = cur_adduct_row[0]
        if mm <= 0:
            return pd.DataFrame()
        else:
            df = pd.DataFrame([[cur_mass_row[0], cur_mass_row[1], mm, mm_low, mm_high, adduct_name]],
                              columns=['ms_name', 'mz', 'observed_mm', 'mm_low', 'mm_high', 'adduct_name'])
            return df

    cur_mass_row = get_mass_range(cur_mass_row, tolerance)

    results_df = []
    for row in adducts_df.values:
        temp = get_cur_mm_adducts(row)
        results_df.append(temp)

    results_df = pd.concat(results_df).reset_index(drop=True)

    return results_df


def main():

    adducts_df = pd.read_csv("./docs/ESI-MS-adducts.csv")
    print("Step 1: Adducts file is loaded!")
    df = pd.read_csv("./Wisecaver/test.txt", index_col=None, header=0)
    print("Step 2: Metabolome file is loaded!")
    df.columns = ['ms_name', 'mz']
    df = df.astype({'mz': float})

    results = []
    tolerance = 30

    with Pool(processes=5) as pool:
        for cur_mass_row in df.values:
            # for i, cur_mass_row in df.iterrows():
            pool.apply_async(get_adduct_data, args=(cur_mass_row, adducts_df, tolerance), callback=results.append)
        pool.close()
        pool.join()
    df = pd.concat(results).reset_index(drop=True)
    print("Step 3: Adducts are added to the metabolome file!")
    temp2 = df.head(n=100)
    temp2.to_csv("table3.csv")

    npdb = NPDB_to_pd("NPDB.sqlite")
    print("Step 4: NPDB is converted to a dataframe!")
    results_df = []
    with Pool(processes=4) as pool:
        for cur_mass_row in df.values:
            pool.apply_async(query_NPDB, args=(cur_mass_row, npdb), callback=results_df.append)
        pool.close()
        pool.join()

    results_df = pd.concat(results_df).reset_index(drop=True)
    temp3 = results_df.head(n=100)
    temp3.to_csv("table5.csv")
    print("Step 5: NPDB is queried with the curated adducts dataframe!")

    # write first set of files
    if not results_df.empty:
        print("Step 6: Writing 1st set of files!")
        res_file = os.path.join("./MassNPDB", 'full.results.csv')
        results_df.to_csv(res_file, sep=',', index=False)
        # ^^ ms_name, mz, observed_mm, NPDB_ID, NPDB_mm, source_id, source_name, InChI, InChI_key, SMILES

        res2_file = os.path.join("./MassNPDB", 'NPDB_entries.txt')
        results_df[['NPDB_ID', 'NPDB_mm', 'SMILES']].drop_duplicates('NPDB_ID').to_csv(res2_file, index=False)
        # ^^ NPDB_ID, NPDB_mm, SMILES


    # Write second set of files
    main_cols = ['ms_name', 'NPDB_ID', 'observed_mm', 'SMILES', 'adduct_name', 'InChI']
    if not results_df.empty:
        print("Step 7: Writing 2nd set of files!")
        # for main results, we only want one entry per ms_name_npdb_id
        res_file = os.path.join("./massNPDB", 'predicted_structures.csv')
        results_df[main_cols].to_csv(res_file, sep=',', index=False)
        # ^^ ms_name, NPDB_ID, NPDB_mm, SMILES, InChI

        res_file = os.path.join("./massNPDB", 'observed_mm.csv')
        results_df[['ms_name', 'mz', 'observed_mm']].drop_duplicates().to_csv(res_file, sep=',', index=False)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    print(end_time - start_time)