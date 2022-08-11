#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import numpy as np
from scipy import stats
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

# project specific module
import gizmos


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('transcripts_file')
    parser.add_argument('metabolites_file')
    parser.add_argument('output_folder')
    parser.add_argument('--cutoff', '-c',
                        default=-.7,
                        required=False,
                        type=float,
                        help='Minimum correlation coefficient; Default: 0.7')
    parser.add_argument('--full_output', '-k',
                        default=False,
                        required=False,
                        action='store_true',
                        help='No cutoff')
    parser.add_argument('--mad_filter', '-d', default=False, required=False, action='store_true', help='Removes arrays with a MAD of 0.')
    parser.add_argument('--remove_zeros', '-z', default=False, required=False, action='store_true', help='Removes arrays with at least one 0.')
    parser.add_argument('--method', '-m', default='spearman', type=str, required=False, choices=['spearman', 'pearson', 'pearsonlog'], help='Default is the spearman correlation method; pearsonlog uses log-transformed arrays, where arrays with atleast a zero are always removed first regardless of otheroptions.')
    parser.add_argument('--annotation', '-a', default='', type=str, required=False, help='Comma-delimited two-columns file with annotations. No header.')
    parser.add_argument('--plot', '-p', default=False, required=False, action='store_true', help='Plots showing the correlation will be created.')
    parser.add_argument('--threads', '-t',
                        default=4,
                        type=int,
                        required=False)
    parser.add_argument('--verbose', '-v',
                        default=False,
                        action='store_true',
                        required=False)
    return parser.parse_args()


# Kumar (2nd Dec 2021)
# applies correlation methods and calculates correlation and
# then merges the results in transcriptomic series as separate columns.
def calc_corr(transcript_series, metabolite_series):
    correlation = 0
    P = 0
    if Options.method == 'pearson':
        correlation, P = stats.pearsonr(transcript_series, metabolite_series)
    elif Options.method == 'pearsonlog':
        correlation, P = stats.pearsonr(np.log10(transcript_series), np.log10(metabolite_series))
    elif Options.method == 'spearman':
        correlation, P = stats.spearmanr(transcript_series, metabolite_series)
    transcript_series['correlation'] = correlation
    transcript_series['P'] = P
    return transcript_series


# Kumar (2nd Dec 2021)
#
def corr_w_transcriptome(metabolite):
    Options = get_args()
    cur_metabolite_transcripts_df = transcripts_df.apply(calc_corr, args=(metabolite,), axis=1)
    if not Options.full_output:
        cutoff_mask = abs(cur_metabolite_transcripts_df.correlation) >= Options.cutoff
        cur_metabolite_transcripts_df = cur_metabolite_transcripts_df[cutoff_mask]

    cur_metabolite_transcripts_df = cur_metabolite_transcripts_df.sort_values(by='correlation', ascending=False)

    results_df = pd.DataFrame({'gene': cur_metabolite_transcripts_df.index, 'metabolite': metabolite.name,
                               'correlation': cur_metabolite_transcripts_df.correlation,
                               'P': cur_metabolite_transcripts_df.P})
    if Options.annotation:
        results_df['gene_annotation'] = annotations_df.annotation[results_df.gene]

    if Options.plot:
        del cur_metabolite_transcripts_df['correlation']
        del cur_metabolite_transcripts_df['P']
        return results_df, metabolite, cur_metabolite_transcripts_df
    else:
        return results_df


def write_and_plot_output(results):
    if Options.plot:
        results_df, metabolite, cur_metabolite_transcripts_df = results
        if not cur_metabolite_transcripts_df.empty:
            # SAVE PLOTS
            if Options.plot:
                # normalizing
                if Options.method == 'pearsonlog':
                    metabolite_n = np.log(metabolite)
                    transcripts_df_n = np.log(cur_metabolite_transcripts_df)
                else:
                    metabolite_n = metabolite.div(metabolite.max())
                    transcripts_df_n = cur_metabolite_transcripts_df.div(cur_metabolite_transcripts_df.max(axis=1),
                                                                         axis=0)

                filename = os.path.join(Options.plots_folder, metabolite.name + '.png')

                cmap = sns.dark_palette("yellow", as_cmap=True)

                # +2 == + white + meta
                fig = plt.figure(figsize=(gizmos.get_heatmap_visual_params(cur_metabolite_transcripts_df.shape[0] + 2,
                                                                           cur_metabolite_transcripts_df.shape[1])))
                gs = gridspec.GridSpec(nrows=cur_metabolite_transcripts_df.shape[0]+2, ncols=1, hspace=0)
                ax1 = fig.add_subplot(gs[0, 0])     # rows, cols
                ax2 = fig.add_subplot(gs[2:, 0])    # rows, cols

                sns.heatmap(data=[metabolite_n], yticklabels=[metabolite_n.name], xticklabels=False,
                            square=True, ax=ax1, cmap=cmap, cbar=False)
                sns.heatmap(data=transcripts_df_n, yticklabels=True, xticklabels=True,
                            square=True, ax=ax2, cmap=cmap, cbar=False)
                ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0)

                # to write annotations, make an invisible plot only with axes on the right.
                # Kumar 04/01/2022
                # to add annotation on the right side of heatmap
                if Options.annotation:
                    n_ticks = len(annotations_df.annotation[transcripts_df_n.index])
                    ticks_labels = annotations_df.annotation[transcripts_df_n.index].iloc[::-1]  # inverted
                    fig = gizmos.plot_annotations(fig, ax2, n_ticks, ticks_labels)

                fig.savefig(filename, bbox_inches='tight')
                plt.close()
    else:
        results_df = results

    # WRITE OUTPUT
    # metabolite, gene, correlation, P, gene_annotation
    with open(Options.output_file, 'a') as f:
        for i, cur_row in results_df.iterrows():
            f.write(cur_row.metabolite + ',')
            f.write(cur_row.gene + ',')
            f.write(str(cur_row.correlation) + ',')
            f.write(str(cur_row.P))
            if Options.annotation:
                f.write(',' + str(cur_row.gene_annotation))
            f.write('\n')
    return


def print_child_proc_error(error_string):
    print('Child process encountered the following error: ' + str(error_string))
    return


# def run_for_multiprocessing():


def main():

    global transcripts_df, annotations_df

    # PARSE INPUT
    Options.output_file = os.path.join(Options.output_folder, 'results.csv')
    Options.log_file = os.path.join(Options.output_folder, 'log.txt')

    # OUTPUT INIT
    gizmos.log_init(Options)

    with open(Options.output_file, 'w') as f:
        if Options.annotation:
            f.write('metabolite,gene,correlation,P,gene_annotation\n')
        else:
            f.write('metabolite,gene,correlation,P\n')

    if Options.plot:
        Options.plots_folder = os.path.join(Options.output_folder, 'plots/')
        if not os.path.exists(Options.plots_folder):
            os.makedirs(Options.plots_folder)

    # LOAD INPUT
    gizmos.print_milestone('Loading input...', Options.verbose)
    transcripts_df = pd.read_csv(Options.transcripts_file, index_col=0)
    metabolites_df = pd.read_csv(Options.metabolites_file, index_col=0)
    if Options.annotation:
        annotations_df = pd.read_csv(Options.annotation, index_col=0, header=None, names=['annotation'])
    else:
        annotations_df = pd.DataFrame({'annotation': []})
    annotations_df.index.name = 'gene'

    for i, row in transcripts_df.iterrows():
        if i not in annotations_df.index:
            annotations_df.annotation[i] = ''

    # MAD FILTER
    if Options.mad_filter:
        transcripts_df = gizmos.apply_MAD_filter(transcripts_df, Options)
        metabolites_df = gizmos.apply_MAD_filter(metabolites_df, Options)

    # REMOVE ZEROES
    if Options.method == 'pearsonlog' or Options.remove_zeros:
        transcripts_df = transcripts_df[(transcripts_df != 0).all(axis=1)]
        metabolites_df = metabolites_df[(metabolites_df != 0).all(axis=1)]

    # SORTING COLUMNS
    gizmos.print_milestone('Preprocessing...', Options.verbose)
    transcripts_labels = set(transcripts_df.columns)
    metabolites_labels = set(metabolites_df.columns)
    common_labels = sorted(list(transcripts_labels.intersection(metabolites_labels)))
    transcripts_df = transcripts_df[common_labels]
    metabolites_df = metabolites_df[common_labels]

    # CORRELATING
    gizmos.print_milestone('Correlating...', Options.verbose)
    if Options.threads == 1:
        for i, cur_metabolite in metabolites_df.iterrows():
            results = corr_w_transcriptome(cur_metabolite)
            write_and_plot_output(results)
    else:
        with Pool(processes=Options.threads) as pool:
            for i, cur_metabolite in metabolites_df.iterrows():
                pool.apply_async(corr_w_transcriptome, args=(cur_metabolite,), callback=write_and_plot_output, error_callback=print_child_proc_error)
            pool.close()
            pool.join()

    return


if __name__ == "__main__":
    # global Options, transcripts_df, annotations_df
    Options = get_args()
    transcripts_df = pd.DataFrame()
    annotations_df = pd.DataFrame()
    main()