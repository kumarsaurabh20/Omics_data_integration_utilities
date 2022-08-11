import io
import os
import sqlite3
import pandas as pd
import numpy as np

#con = sqlite3.connect("mvc.db")
con = sqlite3.connect("NPDB_merge.sqlite")
cur = con.cursor()

#sql = "SELECT uniprot_clusters.uniprot_id, rules.reaction_id " \
#      "FROM uniprot_clusters " \
#      "INNER JOIN smiles_clusters ON smiles_clusters.cluster_id = uniprot_clusters.cluster_id " \
#      "INNER JOIN rules ON rules.smiles_id = smiles_clusters.smiles_id;"

sql = "SELECT * FROM correlations_HB;"
cursor = cur.execute(sql)
temp = []
for each in cursor:
    print(each)
    temp.append(each)

con.close()
uniprot_ids = pd.DataFrame(temp)
uniprot_ids.to_csv("HB_correlations.tsv", sep="\t")

