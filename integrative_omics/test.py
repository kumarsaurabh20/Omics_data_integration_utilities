import pandas as pd
import numpy as np
import gizmos
# save numpy array as csv file
from numpy import asarray
from numpy import savetxt


annotations_df = pd.read_csv("gene_annotation.csv", index_col=None)
annotations_df = annotations_df.rename(columns={annotations_df.columns[0]: 'gene', annotations_df.columns[1]: 'enzyme_pfams'})
enzyme_pfams_list = annotations_df.enzyme_pfams.apply(gizmos.pd_to_list, separator=';')
# expand gene annotations so there's one pfam per line (but we keep the "pfam" annotation that have them all)
lens = [len(item) for item in enzyme_pfams_list]
new_df = pd.DataFrame({'gene': np.repeat(annotations_df.gene, lens), 'pfam_rule': np.concatenate(enzyme_pfams_list)})
annotations_df = pd.merge(annotations_df, new_df, how='outer')



pfam_rules_df = pd.read_csv("pfam_RR_annotation_file.csv", index_col=None)
pfam_rules_df = pfam_rules_df.rename(columns={pfam_rules_df.columns[1]: 'uniprot_id', pfam_rules_df.columns[2]: 'uniprot_enzyme_pfams_acc'})

pfam_rules_df['reaction_id'] = pfam_rules_df.reaction_id.astype('str')



rule = "loose"
        # filter type of anotations (strict, medium, loose)
if rule == 'strict':
    pfam_rules_df = pfam_rules_df[pfam_rules_df.experimentally_validated]
elif rule == 'medium':
    pfam_rules_df = pfam_rules_df[pfam_rules_df.experimentally_validated | pfam_rules_df.Pfam_ec_prediction]
else:  # loose
    pass  # they are all there

# convert pfam_acc to pfam
pfam_dict = pd.read_csv("pfams_dict.csv", index_col=None)
pfam_dict.index = pfam_dict.Acc.apply(lambda x: x.split('.')[0])  # Acc,Name,Desc

uniprot_enzyme_pfams_acc_list = pfam_rules_df.uniprot_enzyme_pfams_acc.apply(gizmos.pd_to_list, separator=';')

pfam_rules_df['uniprot_enzyme_pfams_list'] = [[k for k in row if k in pfam_dict.index] for row in uniprot_enzyme_pfams_acc_list]
pfam_rules_df['uniprot_enzyme_pfams'] = pfam_rules_df.uniprot_enzyme_pfams_list.apply(';'.join)


# Expand df so there is only one pfam per row.
lens = [len(item) for item in pfam_rules_df.uniprot_enzyme_pfams_list]
pfams_flat = [n for sub in pfam_rules_df.uniprot_enzyme_pfams_list for n in sub]
new_df = pd.DataFrame({'uniprot_id': np.repeat(pfam_rules_df.uniprot_id, lens), 'pfam_rule': pfams_flat}).drop_duplicates()
pfam_rules_df = pd.merge(pfam_rules_df, new_df, how='outer')  # on uniprot_id.

del pfams_flat, uniprot_enzyme_pfams_acc_list, new_df
del pfam_rules_df['uniprot_enzyme_pfams_list'], pfam_rules_df['uniprot_enzyme_pfams_acc']

merged_df = pd.merge(annotations_df, pfam_rules_df, how='inner', on="pfam_rule")
# ++ pfam_rule, uniprot_enzyme_pfams >>> -- enzyme_pfams_acc

correlation_df = pd.read_csv("MEJA_correlations.csv", index_col=None)
correlation_df = correlation_df.rename(columns={correlation_df.columns[0]: 'ms_name', correlation_df.columns[1]: 'gene', correlation_df.columns[2]: 'correlation', correlation_df.columns[3]: 'P'})

merged_df = pd.merge(merged_df, correlation_df, how='inner', on='gene')

enzyme_df = merged_df[['reaction_id', 'ms_name', 'gene']].drop_duplicates()

print(enzyme_df)
print(enzyme_df.columns)

transitions_df = pd.read_csv("./test_db/format_database/MassTransitions.csv", index_col=None)
transitions_df['reaction_id'] = transitions_df.reaction_id.astype('str')
transitions_df['mass_transition_round_abs'] = transitions_df.mass_transition_round.apply(abs)
transitions_df = transitions_df[transitions_df.mass_transition_round_abs > 0]
del transitions_df['mass_transition_round_abs']

enzyme_transition_df = pd.merge(enzyme_df, transitions_df, on="reaction_id")
print(enzyme_transition_df)
print(enzyme_transition_df.columns)

enzyme_transition_df = enzyme_transition_df.drop(columns=['gene', 'ms_name']).drop_duplicates()
print(enzyme_transition_df)
print(enzyme_transition_df.columns)

unique_transitions = enzyme_transition_df.mass_transition.unique()
print(unique_transitions)
print(len(unique_transitions))

masses_df = pd.read_csv("./Wisecaver/Wisecaver/merged_metabolome_MeJA.csv", index_col=None, header=0)
masses_df = masses_df.rename(columns={masses_df.columns[0]: 'ms_name', masses_df.columns[1]: 'mz', masses_df.columns[2]: 'mm'}).set_index('ms_name')
print(masses_df)
print(len(masses_df))

path_transitions = np.vstack(masses_df.mm.to_numpy()) + unique_transitions
# define data
#savetxt('path_transitions.csv', path_transitions, delimiter=',')