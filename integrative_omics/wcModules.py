#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import os
import sys
import argparse
import subprocess as sp
from multiprocessing import Pool
from itertools import combinations

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from scipy import stats

import gizmos


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file')
    parser.add_argument('output_folder')
    parser.add_argument('--input_type', '-i',
                        default='transcripts',
                        required=False,
                        type=str,
                        choices=['transcripts', 'correlation', 'MR', 'edges', 'clusterone'],
                        help='Transcripts in RPKM csv and format: rows=genes, columns=conditions'
                             '\nCorrelation in three columns csv with header and format: gene1,gene,PCC'
                             '\nMR in three columns csv with header and format: gene1,gene2,MR'
                             '\nEdges in three columns csv with header and format: gene1,gene2,weight'
                             '\nClusterone output in csv format')
    parser.add_argument('--clusterone_folder', '-1',
                        default='./',
                        required=False,
                        type=str,
                        help='Specify the location of the clusterone jar or an alias.'
                             '\nDefault is at the same folder as the input.')
    parser.add_argument('--edge_weight_cutoff', '-w',
                        default=0.01,
                        required=False,
                        type=float,
                        help='Minimum MR-derived edge weight. Default: 0.01')
    parser.add_argument('--corr_cutoff', '-c',
                        default=0.3,
                        required=False,
                        type=float,
                        help='Minimum correlation coefficient. Default: 0.3')
    parser.add_argument('--mad_filter', '-d',
                        default=False,
                        required=False,
                        action='store_true',
                        help='Removes arrays with a MAD of 0.')
    parser.add_argument('--remove_zeros', '-z',
                        default=False,
                        required=False,
                        action='store_true',
                        help='Removes arrays with at least one 0.')
    parser.add_argument('--method', '-m',
                        default='pearson',
                        type=str,
                        required=False,
                        choices=['pearson', 'pearsonlog'],
                        help='Default: pearson. pearsonlog uses log-transformed arrays, for which arrays with any 0 '
                             'are always removed first, regardless of other options.')
    parser.add_argument('--annotations_file', '-a',
                        default='',
                        type=str,
                        required=False,
                        help='Comma-delimited two-columns file with annotations. No header.')
    parser.add_argument('--expression', '-x',
                        default='',
                        type=str,
                        required=False,
                        help='Only use this option for plotting when initial input is not transcripts. '
                             'Transcripts in RPKM csv and format: rows=genes, columns=conditions')
    parser.add_argument('--threads', '-t',
                        default=4,
                        type=int,
                        required=False)
    parser.add_argument('--plots', '-p',
                        default=False,
                        action='store_true',
                        required=False,
                        help='Will produce plots for resulting wc modules.')
    parser.add_argument('--verbose', '-v',
                        default=False,
                        action='store_true',
                        required=False)
    return parser.parse_args()


def calc_MR(rank1, rank2):
    MR = np.sqrt(rank1 * rank2)
    return MR


def calc_corr(transcript_series1, transcript_series2):
    method = 'pearson'
    correlation = 0
    P = 0
    if method == 'pearson':
        correlation, P = stats.pearsonr(transcript_series1, transcript_series2)
        print(correlation)
    else:   # pearsonlog
        correlation, P = stats.pearsonr(np.log10(transcript_series1), np.log10(transcript_series2))
        print(correlation)

    if correlation >= 0.3:
        result_string = ','.join([transcript_series1.name, transcript_series2.name, str(correlation)])
    else:
        result_string = ''

    print(result_string)
    return result_string


def get_annotation_members(cur_members, annotations_series):
    annotations = annotations_series[cur_members]
    annotations = [n.split(';') for n in annotations]
    annotations = {n for sub in annotations for n in sub}
    if '' in annotations:
        annotations.remove('')
    return annotations


def write_correlation(result_string):
    if result_string:
        with open(Options.correlations_file, 'a') as f:
            f.write(result_string)
            f.write('\n')
    return


def get_correlations_df_from_transcripts():
    # LOAD INPUT
    gizmos.print_milestone('Loading input...', Options.verbose)
    transcripts_df = pd.read_csv(Options.input_file, index_col=0)
    print(transcripts_df.shape)
    transcripts_df.index.name = 'gene'


    # MAD FILTER
    if Options.mad_filter:
        transcripts_df = gizmos.apply_MAD_filter(transcripts_df, Options)

    # REMOVE ZEROES
    if Options.method == 'pearsonlog' or Options.remove_zeros:
        transcripts_df = transcripts_df[(transcripts_df != 0).all(axis=1)]

    print(transcripts_df.index)

    # CORRELATING
    gizmos.print_milestone('Correlating...', Options.verbose)
    with Pool(processes=Options.threads) as pool:
        for transcript1, transcript2 in combinations(transcripts_df.index, 2):
            pool.apply_async(calc_corr, args=(transcripts_df.loc[transcript1], transcripts_df.loc[transcript2]), callback=write_correlation)
        pool.close()
        pool.join()

    return


def get_MR_from_correlations(correlations_df, all_genes):
    # GET RANKS
    gizmos.print_milestone('Getting ranks...', Options.verbose)
    ranks_dict = {}
    for cur_gene in all_genes:
        df1 = correlations_df[correlations_df.gene1 == cur_gene]
        df1 = df1[['gene2', 'correlation']]
        df1.columns = ['other_gene', 'correlation']

        df2 = correlations_df[correlations_df.gene2 == cur_gene]
        df2 = df2[['gene1', 'correlation']]
        df2.columns = ['other_gene', 'correlation']

        dff = pd.concat([df1, df2]).sort_values(by='correlation', ascending=False)
        dff.index = np.arange(1, len(dff) + 1)  # this makes index = rank(cur_gene -> other_gene)

        ranks_dict[cur_gene] = {}
        for rank, cur_row in dff.iterrows():
            ranks_dict[cur_gene][cur_row.other_gene] = rank

    # GET MUTUAL RANKS
    gizmos.print_milestone('Getting mutual ranks...', Options.verbose)
    edges_df = []
    gene_combinations = {n for n in combinations(all_genes, 2)}
    for gene1, gene2 in gene_combinations:
        if ((gene1 in ranks_dict and
             gene2 in ranks_dict and
             gene2 in ranks_dict[gene1] and
             gene1 in ranks_dict[gene2])):
            rank1 = ranks_dict[gene1][gene2]
            rank2 = ranks_dict[gene2][gene1]
            edges_df.append([gene1, gene2, calc_MR(rank1, rank2)])
    edges_df = pd.DataFrame.from_records(edges_df, columns=['gene1', 'gene2', 'MR'])

    # OUTPUT
    fname = os.path.join(Options.output_folder, 'MR.csv')
    edges_df.to_csv(fname, index=False)
    return edges_df


def get_weight_from_MR(edges_df):
    # CONVERT TO EDGE WEIGHT
    edges_df['weight'] = np.exp(-(edges_df.MR - 1) / 25)  # other options: /5 /10 /25 /50 /100
    edges_df = edges_df[edges_df.weight >= Options.edge_weight_cutoff].reset_index(drop=True)
    del edges_df['MR']

    gizmos.print_milestone('Writing edges...', Options.verbose)
    Options.edges_file = os.path.join(Options.output_folder, 'edges.txt')
    edges_df.to_csv(Options.edges_file, index=False, header=False, sep=' ')
    return


def run_clusterone():
    # CLUSTERONE
    gizmos.print_milestone('Using ClusterOne...', Options.verbose)
    print('\n')
    with open(Options.c1_file, 'w') as f:
        # java -jar cluster_one-1.0.jar input -F csv
        sp.run(['java', '-jar', Options.clusterone, Options.edges_file, '-F', 'csv'], stdout=f)
    print('\n')
    return


def use_pfam_rules(cur_cluster, pfam_rules_df, annotations_series):
    bio_df = pfam_rules_df.loc[cur_cluster.bio_pfams]
    core_mask = bio_df.enzyme_type == 'Core'
    tailor_mask = bio_df.enzyme_type == 'Other'
    cores = set(bio_df[core_mask].index)
    tailors = set(bio_df[tailor_mask].index)

    if cores:
        core_genes = set()
        for gene in cur_cluster.Members:
            if set(annotations_series[gene].split(';')).intersection(cores):
                core_genes.add(gene)
        cur_cluster['has_core'] = True
        cur_cluster['metabolite_type'] = set(pfam_rules_df.metabolite_type.loc[cores])
        cur_cluster['core_genes'] = core_genes
    else:
        cur_cluster['has_core'] = False
        cur_cluster['metabolite_type'] = set()
        cur_cluster['core_genes'] = set()

    if tailors:
        tailor_genes = set()
        for gene in cur_cluster.Members:
            if set(annotations_series[gene].split(';')).intersection(tailors):
                tailor_genes.add(gene)
        cur_cluster['has_tailoring'] = True
        cur_cluster['tailoring_genes'] = tailor_genes
    else:
        cur_cluster['has_tailoring'] = False
        cur_cluster['tailoring_genes'] = set()

    return cur_cluster


def plot_output(cur_module, full_transcripts_df, annotations_series):
    # cur_module = [Cluster], Size, Members, bio_pfams, metabolite_type, core_genes, tailoring_genes

    genes = gizmos.pd_to_set(cur_module.Members, ';')
    transcripts_df = full_transcripts_df.loc[genes]

    # normalizing
    if Options.method == 'pearsonlog':
        transcripts_df_n = np.log(transcripts_df)
    else:
        transcripts_df_n = transcripts_df.div(transcripts_df.max(axis=1), axis=0)

    # splitting
    core_genes = gizmos.pd_to_set(cur_module.core_genes, ';')
    tailoring_genes = gizmos.pd_to_set(cur_module.tailoring_genes, ';')
    genes.difference_update(core_genes)
    genes.difference_update(tailoring_genes)

    transcripts_core = transcripts_df_n.loc[core_genes]
    transcripts_tailoring = transcripts_df_n.loc[tailoring_genes]
    transcripts_rest = transcripts_df_n.loc[genes]

    ax_core_end = len(transcripts_core)
    ax_tailoring_start = ax_core_end + 1
    ax_tailoring_end = ax_tailoring_start + len(transcripts_tailoring)
    ax_rest_start = ax_tailoring_end + 1

    fname = os.path.join(Options.plots_folder, 'cluster_' + str(cur_module.name) + '.png')

    cmap = sns.dark_palette("yellow", as_cmap=True)

    fig = plt.figure(figsize=(gizmos.get_heatmap_visual_params(transcripts_df.shape[0] + 2, # + 2 whites
                                                               transcripts_df.shape[1])))
    gs = gridspec.GridSpec(nrows=len(transcripts_df_n)+2, ncols=1, hspace=0)

    ax_core = fig.add_subplot(gs[:ax_core_end, 0])     # rows, cols
    ax_tailoring = fig.add_subplot(gs[ax_tailoring_start:ax_tailoring_end, 0])     # rows, cols
    ax_rest = fig.add_subplot(gs[ax_rest_start:, 0])     # rows, cols

    sns.heatmap(data=transcripts_core, yticklabels=True, xticklabels=False,
                square=True, ax=ax_core, cmap=cmap, cbar=False)
    sns.heatmap(data=transcripts_tailoring, yticklabels=True, xticklabels=False,
                square=True, ax=ax_tailoring, cmap=cmap, cbar=False)
    sns.heatmap(data=transcripts_rest, yticklabels=True, xticklabels=True,
                square=True, ax=ax_rest, cmap=cmap, cbar=False)

    ax_core.set_yticklabels(ax_core.get_yticklabels(), rotation=0)
    ax_tailoring.set_yticklabels(ax_tailoring.get_yticklabels(), rotation=0)
    ax_rest.set_yticklabels(ax_rest.get_yticklabels(), rotation=0)

    # to write annotations, make an invisible plot only with axes on the right, for each plot
    for cur_ax, cur_transcripts in [[ax_core, transcripts_core],
                                    [ax_tailoring, transcripts_tailoring], [ax_rest, transcripts_rest]]:
        n_ticks = len(annotations_series[cur_transcripts.index])
        ticks_labels = annotations_series[cur_transcripts.index].iloc[::-1]  # inverted
        fig = gizmos.plot_annotations(fig, cur_ax, n_ticks, ticks_labels)

    fig.savefig(fname, bbox_inches='tight')
    plt.close()
    return


#################
# MAIN PIPELINE #
#################
def main():
    # PARSE INPUT
    Options.c1_file = os.path.join(Options.output_folder, 'clusterone_output.csv')
    Options.log_file = os.path.join(Options.output_folder, 'log.txt')

    # CHECKING INPUTS EXIST
    # TODO make actual errors
    if not os.path.exists(Options.input_file):
        gizmos.print_milestone('Error: could not locate input_file', Options.verbose)
    Options.clusterone = os.path.join(Options.clusterone_folder, 'cluster_one-1.0.jar')
    if not os.path.exists(Options.clusterone):
        gizmos.print_milestone('Warning: could not locate cluster_one-1.0.jar', Options.verbose)
    if Options.expression:
        if not os.path.exists(Options.expression):
            gizmos.print_milestone('Error: could not locate the expression file.', Options.verbose)

    # OUTPUT INIT
    gizmos.log_init(Options)

    # ### INPUT = TRANSCRIPTS
    if Options.input_type == 'transcripts':
        # PREPARE OUTPUT AND VARIABLES
        Options.expression = Options.input_file
        Options.correlations_file = os.path.join(Options.output_folder, 'correlations.csv')
        with open(Options.correlations_file, 'w') as f:
            f.write('gene1,gene2,correlation')
            f.write('\n')
        # PIPELINE
        get_correlations_df_from_transcripts()
        correlations_df = pd.read_csv(Options.correlations_file, index_col=None)
        correlations_df.columns = ['gene1', 'gene2', 'correlation']
        all_genes = set(correlations_df.gene1).union(correlations_df.gene2)
        edges_df = get_MR_from_correlations(correlations_df, all_genes)
        get_weight_from_MR(edges_df)
        run_clusterone()

    # ### INPUT = CORRELATION
    elif Options.input_type == 'correlation':
        # LOAD INPUT
        gizmos.print_milestone('Loading input...', Options.verbose)
        correlations_df = pd.read_csv(Options.input_file, index_col=None)
        correlations_df.columns = ['gene1', 'gene2', 'correlation']
        # CORRELATION CUTOFF
        mask = correlations_df.correlation >= Options.corr_cutoff
        correlations_df = correlations_df[mask]
        all_genes = set(correlations_df.gene1).union(correlations_df.gene2)
        # PIPELINE
        edges_df = get_MR_from_correlations(correlations_df, all_genes)
        get_weight_from_MR(edges_df)
        run_clusterone()

    # ### INPUT = MR
    elif Options.input_type == 'MR':
        # LOAD INPUT
        gizmos.print_milestone('Loading input...', Options.verbose)
        edges_df = pd.read_csv(Options.input_file, index_col=None)
        edges_df.columns = ['gene1', 'gene2', 'MR']
        # PIPELINE
        get_weight_from_MR(edges_df)
        run_clusterone()

    # ### INPUT = EDGES
    elif Options.input_type == 'edges':
        Options.edges_file = Options.input_file
        run_clusterone()

    # ### INPUT = CLUSTERONE
    elif Options.input_type == 'clusterone':
        Options.c1_file = Options.input_file

    # ### APPLY PFAM RULES
    if Options.annotations_file:
        # READ INPUT
        gizmos.print_milestone('Reading ClusterOne results and annotations...', Options.verbose)
        Options.pfam_rules_file = os.path.join(sys.path[0], 'psmash_rules.csv')
        pfam_rules_df = pd.read_csv(Options.pfam_rules_file, index_col=0)
        all_bio_pfams = set(pfam_rules_df.index)

        annotations_series = pd.read_csv(Options.annotations_file, index_col=0, header=None, names=['annotation'], squeeze=True)

        c1_df = pd.read_csv(Options.c1_file, index_col=0)
        pmask = c1_df['P-value'] <= 0.1
        c1_df = c1_df[pmask]
        c1_df['Members'] = c1_df.Members.apply(gizmos.pd_to_set, args=(' ',))
        all_genes = set.union(*c1_df.Members)

        # ANNOTATE MISSING GENES IN ORIGINAL ANNOTATIONS
        for gene in all_genes:
            if gene not in annotations_series:
                annotations_series[gene] = ''

        # GET CLUSTER PFAMS
        c1_df['pfams'] = c1_df.Members.apply(get_annotation_members, args=(annotations_series,))
        c1_df['bio_pfams'] = c1_df.pfams.apply(all_bio_pfams.intersection)

        # APPLY PFAM RULES
        c1_df = c1_df.apply(use_pfam_rules, args=(pfam_rules_df, annotations_series,), axis=1)

        # OUTPUT
        gizmos.print_milestone('Writing output...', Options.verbose)
        c1_df['Members'] = c1_df.Members.apply(gizmos.set_to_string, args=(';',))
        c1_df['pfams'] = c1_df.pfams.apply(gizmos.set_to_string, args=(';',))
        c1_df['bio_pfams'] = c1_df.bio_pfams.apply(gizmos.set_to_string, args=(';',))
        c1_df['metabolite_type'] = c1_df.metabolite_type.apply(gizmos.set_to_string, args=(';',))
        c1_df['core_genes'] = c1_df.core_genes.apply(gizmos.set_to_string, args=(';',))
        c1_df['tailoring_genes'] = c1_df.tailoring_genes.apply(gizmos.set_to_string, args=(';',))
        fname = os.path.join(Options.output_folder, 'clusterone.annotated.csv')
        c1_df.to_csv(fname, index=True)

        fname = os.path.join(Options.output_folder, 'wcmodules.csv')
        c1_df = c1_df[c1_df.has_core & c1_df.has_tailoring]
        c1_df = c1_df[['Size', 'Members', 'bio_pfams', 'metabolite_type', 'core_genes', 'tailoring_genes']]
        c1_df.to_csv(fname, index=True)

        if Options.plots:
            # CHECKING INPUTS
            if not Options.expression:
                gizmos.print_milestone('Error: plotting flag was added, but not an expression file. Cannot plot.',
                                       Options.verbose)
            else:
                transcripts_df = pd.read_csv(Options.expression, index_col=0)
                transcripts_df.index.name = 'gene'

                # PLOTS FOLDER
                Options.plots_folder = os.path.join(Options.output_folder, 'plots/')
                if not os.path.exists(Options.plots_folder):
                    os.makedirs(Options.plots_folder)

                # PLOTTING
                gizmos.print_milestone('Plotting...', Options.verbose)
                c1_df.apply(plot_output, args=(transcripts_df, annotations_series,), axis=1)

    return


if __name__ == "__main__":
    Options = get_args()
    main()
