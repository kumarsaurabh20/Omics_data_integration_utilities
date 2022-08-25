import io
import os
import sqlite3
import pandas as pd
import numpy as np

#con = sqlite3.connect("mvc.db")
con = sqlite3.connect("/Users/singh018/Documents/NPDB_merge.sqlite")
cur = con.cursor()

#sql = "SELECT uniprot_clusters.uniprot_id, rules.reaction_id " \
#      "FROM uniprot_clusters " \
#      "INNER JOIN smiles_clusters ON smiles_clusters.cluster_id = uniprot_clusters.cluster_id " \
#      "INNER JOIN rules ON rules.smiles_id = smiles_clusters.smiles_id;"

sql = "SELECT ms_name,mz, monoisotopic_mass, smile, inchi_key2 FROM merged_data;"
cursor = cur.execute(sql)
temp = []
for each in cursor:
    print(each)
    temp.append(each)

con.close()
uniprot_ids = pd.DataFrame(temp)
uniprot_ids.columns = ["ms_name", "mz", "mm", "smile", "inchi_key2"]
print("Original count of the dataframe without filtering {}".format(uniprot_ids.shape))
final = uniprot_ids.drop_duplicates(subset=['inchi_key2'], keep="first")
print("Shape of the dataframe with inchi_key2 as filtering criteria {}".format(final.shape))
final2 = uniprot_ids.drop_duplicates(subset=['smile'], keep="first")
print("Shape of the dataframe with smile as filtering criteria {}".format(final2.shape))
new_masses = final2[['ms_name', 'mz', 'mm']]
new_masses.to_csv("smile_queryMassNPDB.csv", sep=",", index=False)

#uniprot_ids.to_csv("queryMassNPDB_hp.csv", sep=",", index=False)
