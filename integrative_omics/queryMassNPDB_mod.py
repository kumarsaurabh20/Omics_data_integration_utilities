#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import time
import numpy as np
import pandas as pd
import sqlite3
#import janitor

import gizmos

def get_args():
    """
    Self-explanatory. Toodles.
    :return:
    """
    parser = argparse.ArgumentParser(prog='queryNPDB',
                                     description='Integrate NPDB/LOTUS and extract structures based on your mass signatures!',
                                     epilog='Contact - kumar.singh@wur.nl')
    parser.add_argument('-add', '--adducts_file', default=True, action='store', required=True, help='ESI-MS_adducts.csv')
    parser.add_argument('-ms', '--mass_signatures_file', default=True, action='store', required=True, help='Mass signatures two column file: {id,mz}')
    parser.add_argument('-npdb', '--npdb_sqlite_file', default=True, action='store', required=True, help='Natural Product db file in sqlite format!')
    parser.add_argument('-dn', '--sqlite_db_name', default=True, action='store', required=True, help='Provide a name for the database!')
    parser.add_argument('-c', '--chunk_size', action='store', default=3, type=int, required=True, help='Divide the mass signatures file in to chunks!')
    parser.add_argument('-p', '--ppm', action='store', default=30, type=int, required=True, help='Divide the input file in to chunks!')
    parser.add_argument('-o', '--output_folder', default=True, action='store', required=True)
    parser.add_argument('-f', '--output_csv_files', default=False, action='store_true', required=False, help='Output integration results in to multiple CSV files!')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', required=False)
    return parser.parse_args()

def split_dataframe(df, chunk_size = 3):
    num_chunks = len(df) // chunk_size
    if len(df) % chunk_size != 0:
        num_chunks += 1
    for i in range(num_chunks):
        yield df[i*chunk_size:(i + 1) * chunk_size]


def NPDB_to_pd(npdb):
    conn = sqlite3.connect(npdb)
    c = conn.cursor()
    query = c.execute("SELECT structure.structure_id,structure.monoisotopic_mass, structure_has_data_source.source_id, structure_has_data_source.source_name, structure.inchi,structure.inchi_key2,structure.smile FROM structure left join structure_has_data_source on structure_has_data_source.structure_id = structure.structure_id")
    cols = [column[0] for column in query.description]
    results = pd.DataFrame.from_records(data=query.fetchall(), columns=cols)
    c.close()

    return results

def checkP(x):
    return x + x*Options.ppm/1000000
def checkN(x):
    return x - x*Options.ppm/1000000

def create_output_files(sqlitedb, sql, columns, filename):
    db_file = os.path.join(sqlitedb)
    con = sqlite3.connect(db_file)
    cur = con.cursor()
    cursor = cur.execute(sql)
    count = 1

    signatures = {}
    for each in cursor:
        #print(each)
        #ms_name, mz, mm = each
        signatures[count] = list(each)
        count += 1
    con.close()

    result = pd.DataFrame(signatures.values(), columns=columns)
    result.to_csv(filename, index=False)


def main():

    # Load the adducts file
    adducts_df = pd.read_csv(Options.adducts_file)
    gizmos.print_milestone('Step 1: Adducts file is loaded!', Options.verbose)

    # Load the mass signature file
    df = pd.read_csv(Options.mass_signatures_file, index_col=None, header=0)
    gizmos.print_milestone('Step 2: Metabolome file is loaded!', Options.verbose)

    # Load the NPDB database
    npdb = NPDB_to_pd(Options.npdb_sqlite_file)

    # Create a folder to write all files related with this script.
    if not os.path.exists(Options.output_folder):
        os.makedirs(Options.output_folder)

    # Create an empty sqlite db
    db_file_path = os.path.join(Options.output_folder, Options.sqlite_db_name)
    conn = sqlite3.connect(db_file_path)
    #
    df.columns = ['ms_name', 'mz']
    df = df.astype({'mz': float})
    chunks = split_dataframe(df, 3)
    counter = 1
    gizmos.print_milestone('Step 3: Processing mass signature file...', Options.verbose)
    for each in chunks:
        gizmos.print_milestone('    Reading chunk {}'.format(counter), Options.verbose)
        chunk_df = pd.DataFrame(each)
        temp = chunk_df.iloc[:, 1]
        ids = chunk_df.iloc[:, 0]

        # transform the mz with +/- ppm
        col1 = temp.transform(checkP)
        col2 = temp.transform(checkN)

        # add +/- ppm to the dataframe
        chunk_df.loc['top'] = col1
        chunk_df.loc['bottom'] = col2

        # make adducts dataframe comparable to mz i.e add all adducts to the individual feature IDs
        new = pd.concat([adducts_df]*len(chunk_df), keys=ids)
        new.reset_index(inplace=True)

        # merge the chunk_df with all all adducts data frame using feature id
        df_outer = pd.merge(chunk_df, new, on="ms_name", how="outer")

        # Create high and low columns
        df_outer['mm'] = df_outer["mz"] #- df_outer["Mass"] / df_outer["Mult"]
        df_outer['mm_high'] = df_outer["mz"] + df_outer["Mass"] / df_outer["Mult"]
        df_outer['mm_low'] = df_outer["mz"] - df_outer["Mass"] / df_outer["Mult"]

        a = npdb.monoisotopic_mass.values
        bh = df_outer.mm_high.values
        bl = df_outer.mm_low.values

        # numpy broadcasting (https://stackoverflow.com/questions/44367672/best-way-to-join-merge-by-range-in-pandas)
        # We look for every instance of a being greater than or equal to bl
        # while at the same time a is less than or equal to bh.
        i, j = np.where((a[:, None] >= bl) & (a[:, None] <= bh))
        ##
        gizmos.print_milestone('    Extracting structures from NPDB...', Options.verbose)

        # we used pandas approach instead of numpy in the generation of final dataframe.
        # This works much faster and consumes less memory.
        # append method linked with concat can get all the NaN values as well from the unmatched rows.
        temp2 = pd.concat([npdb.loc[i, :].reset_index(drop=True), df_outer.loc[j, :].reset_index(drop=True)], axis=1).append(npdb[~np.in1d(np.arange(len(npdb)), np.unique(i))],ignore_index=True, sort=False)

        # drop the NaN rows and get rid of the redundancy by droping duplicates based on inchi key.
        final = temp2[temp2['mm_high'].notnull()].drop_duplicates(subset=['inchi_key2'], keep="first")

        # Drop the unwanted columns from the final dataframe
        final.drop(['level_1', 'source_id', 'Ion_mode', 'Ion_mass', 'Charge', 'Mult', 'Mass', 'mm', 'mm_low', 'mm_high'], axis=1, inplace=True)

        # Push the final dataframe in a SQLite database, with append on as we are looping over chunks.
        gizmos.print_milestone('    Writing the intersection to the database...', Options.verbose)
        final.to_sql("merged_data", conn, if_exists = "append" ,index=False)

        counter += 1
    conn.close()

    # Generate NPDB formatted files
    if Options.output_csv_files:
        gizmos.print_milestone('    Writing CSV files...', Options.verbose)
        res_file_path = os.path.join(Options.output_folder, 'full_results.csv')
        sql_full = "SELECT * FROM merged_data;"
        create_output_files(db_file_path, sql_full, ['NPDB_ID', 'mm', 'adduct_name', 'source_name', 'INCHI', 'INCHI_key2', 'SMILES', 'ms_name', 'mz'], res_file_path)

        res_file_path = os.path.join(Options.output_folder, 'NPDB_entries.csv')
        sql_npdb_entries = "SELECT structure_id, monoisotopic_mass, smile FROM merged_data;"
        create_output_files(db_file_path, sql_npdb_entries, ['NPDB_ID','NPDB_mm','SMILES'], res_file_path)

        res_file_path = os.path.join(Options.output_folder, 'observed_mm.csv')
        sql_mm = "SELECT ms_name, mz, monoisotopic_mass FROM merged_data;"
        create_output_files(db_file_path, sql_mm, ['ms_name', 'mz', 'mm'], res_file_path)

        res_file_path = os.path.join(Options.output_folder, 'predicted_structures.csv')
        sql_structures = "SELECT ms_name, structure_id, monoisotopic_mass, smile, Ion_name, inchi FROM merged_data;"
        create_output_files(db_file_path, sql_structures, ['ms_name', 'NPDB_ID', 'mm', 'SMILE', 'adduct_name', 'INCHI'], res_file_path)

if __name__ == "__main__":
    Options = get_args()
    start_time = time.time()
    main()
    end_time = time.time()
    print(end_time - start_time)

#def create_connection(db_file):
#    """ create a database connection to a SQLite database """
#    conn = None
#    try:
#        conn = sqlite3.connect(db_file)
#        print(sqlite3.version)
#    except Error as e:
#        print(e)
#    finally:
#        if conn:
#            conn.close()
# if __name__ == '__main__':
#   create_connection(r"C:\sqlite\db\pythonsqlite.db")

#npdb.where(npdb['monoisotopic_mass']>=df_outer['mm_low']).where(npdb['A']<df_outer['mm_high']).dropna()
        #can only compare identically labelled series objects

        #final = npdb.conditional_join(df_outer,('monoisotopic_mass', 'mm_low', '>='),('monoisotopic_mass', 'mm_high', '<='),how = 'left')
        #print(final[final['mm_high'].notnull()])
        #janitor based conditional join statement
        #print(df_outer)

        #SQLlite based solution
        #conn = sqlite3.connect(":memory:")
        #npdb.to_sql("npdb", conn, index=False)
        #df_outer.to_sql("df_outer", conn, index=False)
        #qry = "SELECT * FROM npdb, df_outer WHERE npdb.monoisotopic_mass >= df_outer.mm_low and npdb.monoisotopic_mass <= df_outer.mm_high"
        #tt = pd.read_sql_query(qry,conn)
        #print(tt)

# SQLlite based solution
        # conn = sqlite3.connect(":memory:")
        # npdb.to_sql("npdb", conn, index=False)
        # df_outer.to_sql("df_outer", conn, index=False)
        # qry = "SELECT * FROM npdb, df_outer WHERE npdb.monoisotopic_mass >= df_outer.mm_low and npdb.monoisotopic_mass <= df_outer.mm_high"
        # tt = pd.read_sql_query(qry,conn)
        # print(tt)

        #final.to_csv(filename, index = False)

#print(df_outer)
        #print(df_outer.columns)
        #print(range)

        #df_outer['temp'] = 1
        #npdb['temp'] = 1
        #final = pd.merge(df_outer, npdb, on="temp", how="outer")
        #print(final.columns)

        #conditions = [final['monoisotopic_mass'].ge(final['mm_low']) & final['monoisotopic_mass'].le(final['mm_high'])]
        #choices = [0]
        #final['filter'] = np.select(conditions, choices, default=0)
        #print(final)
        #print(final.columns)

#def split_dataframe(df, chunk_size = 2):
#    chunks = list()
#    num_chunks = math.ceil(len(df)) / chunk_size
#    for i in range(int(num_chunks)):
#        chunks.append(df[i*chunk_size:(i+1)*chunk_size])
#    return chunks