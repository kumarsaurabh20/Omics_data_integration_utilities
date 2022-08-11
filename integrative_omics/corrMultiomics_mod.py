#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import numpy as np
from scipy import stats
from multiprocessing import Pool, Manager
import sqlite3
from multiprocessing.shared_memory import SharedMemory
from multiprocessing.managers import SharedMemoryManager
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

# project specific module
#import gizmos

def calc_corr(transcript_series, metabolite_series):
    method = "pearson"
    correlation = 0
    P = 0
    if method == 'pearson':
        correlation, P = stats.pearsonr(transcript_series, metabolite_series)
    elif method == 'pearsonlog':
        correlation, P = stats.pearsonr(np.log10(transcript_series), np.log10(metabolite_series))
    elif method == 'spearman':
        correlation, P = stats.spearmanr(transcript_series, metabolite_series)
    transcript_series['correlation'] = correlation
    transcript_series['P'] = P
    return transcript_series


def write_and_plot_output(results):
    results_df = results
    # WRITE OUTPUT
    # metabolite, gene, correlation, P, gene_annotation
    output_file = "results.csv"

    with open(output_file, 'a') as f:
        for i, cur_row in results_df.iterrows():
            f.write(cur_row.metabolite + ',')
            f.write(cur_row.gene + ',')
            f.write(str(cur_row.correlation) + ',')
            f.write(str(cur_row.P) + "\n")
    return


def corr_w_transcriptome(metabolite, transcripts_df):
    cur_metabolite_transcripts_df = transcripts_df.apply(calc_corr, args=(metabolite,), axis=1)
    # print(cur_metabolite_transcripts_df)
    cutoff_mask = abs(cur_metabolite_transcripts_df.correlation) >= 0.7
    cur_metabolite_transcripts_df = cur_metabolite_transcripts_df[cutoff_mask]
    cur_metabolite_transcripts_df = cur_metabolite_transcripts_df.sort_values(by='correlation', ascending=False)
    # print(cur_metabolite_transcripts_df)
    results_df = pd.DataFrame({'metabolite': metabolite.name, 'gene': cur_metabolite_transcripts_df.index, 'correlation': cur_metabolite_transcripts_df.correlation, 'P': cur_metabolite_transcripts_df.P})
    print(results_df)
    return results_df


def print_child_proc_error(error_string):
    print('Child process encountered the following error: ' + str(error_string))
    return

def calc_MAD(x):
    """
    Returns the median absolute deviation via the following equation:
    MAD = median( | Xi - median(X) |  )
    :param x:
    :return:
    """
    med = np.median(x)
    x2 = np.absolute(x-med)
    MAD = np.median(x2)
    return MAD


def apply_MAD_filter(df):
    df['MAD'] = df.apply(calc_MAD, axis=1)
    mask = df.MAD > 0
    df_filtered = df[mask]
    del df_filtered['MAD']
    return df_filtered


if __name__ == "__main__":
    # global Options, transcripts_df, annotations_df
    pool = Pool(processes=8, maxtasksperchild=10)
    # smm = SharedMemoryManager()
    transcripts_df = pd.read_csv("/Users/singh018/Documents/Meantools_v1/Wisecaver/Trans_quant_HB.csv", index_col=0)
    metabolites_df = pd.read_csv("/Users/singh018/Documents/Meantools_v1/Wisecaver/Metabo_quant_HB.csv", index_col=0)
    print("Step 1")

    transcripts_df = apply_MAD_filter(transcripts_df)
    metabolites_df = apply_MAD_filter(metabolites_df)
    print("Step 2")

    transcripts_df = transcripts_df[(transcripts_df != 0).all(axis=1)]
    metabolites_df = metabolites_df[(metabolites_df != 0).all(axis=1)]
    print("Step 3")

    transcripts_labels = set(transcripts_df.columns)
    metabolites_labels = set(metabolites_df.columns)
    common_labels = sorted(list(transcripts_labels.intersection(metabolites_labels)))
    transcripts_df = transcripts_df[common_labels]
    metabolites_df = metabolites_df[common_labels]
    print("Step 4")

    # Kumar 2nd Jan 2022
    # Manager class to make all process share the global variables.
    mgr = Manager()
    ns = mgr.Namespace()
    ns.df = transcripts_df
    results = []

    #with Pool(processes=8, maxtasksperchild=10) as pool:
    for i, cur_metabolite in metabolites_df.iterrows():
        pool.apply_async(corr_w_transcriptome, args=(cur_metabolite, ns.df), error_callback=print_child_proc_error, callback=results.append)
    pool.close()
    pool.join()

    results_df = pd.concat(results).reset_index(drop=True)
    conn = sqlite3.connect("NPDB_merge.sqlite")
    results_df.to_sql("correlations_HB", conn, if_exists="append", index=False)
    conn.close()
