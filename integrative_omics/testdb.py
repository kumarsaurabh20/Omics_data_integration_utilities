import io
import os
import sqlite3
from Bio.KEGG import REST
from Bio.KEGG import Enzyme
import pandas as pd
import numpy as np
from prody import *
import urllib3

################ Uniprot - Pfam mapping ################
# urllib3.disable_warnings()
# query = ["P40925", "P40926", "O43175", "Q9UM73", "P97793"]
# map = {}
#for id in query:
#   temp = searchPfam(id)
#   temp2 = []
#   for each in temp.keys():
#      temp2.append(each)
#   map[id] = temp2
#print(map)
################ Uniprot - Pfam mapping ################

# Some code to return a Pandas dataframe, given tabular text
def rest_to_df(result):
    temp = pd.read_table(io.StringIO(result), header=None)
    temp[0] = temp[0].str.replace('ko:', '')
    temp[1] = temp[1].str.replace('rn:', '')
    return temp

con = sqlite3.connect("mvc.db")
cur = con.cursor()

# SELECT * FROM chemical_species WHERE kegg IS NOT NULL AND metacyc IS NOT NULL AND reactome IS NOT NULL and chebi IS NOT NULL;
# sql1 = "SELECT name FROM sqlite_schema WHERE type = 'table' AND name NOT LIKE 'sqlite_%'"
# sql2 = "SELECT * FROM '{}'"
# sql = "SELECT kegg FROM reactions WHERE kegg IS NOT NULL AND metacyc IS NOT NULL AND reactome IS NOT NULL and rhea IS NOT NULL;"
sql = "SELECT kegg FROM reactions WHERE kegg IS NOT NULL;"
cursor = cur.execute(sql)
counter = 0
kegg_reaction_ids = []
for each in cursor:
    kegg_reaction_ids.append(''.join(list(each)))
# print(kegg_reaction_ids)
con.close()

#### Query KEGG ####
results = REST.kegg_link("reaction", "ko").read()
reaction_ko_df = rest_to_df(results)
reaction_ko_df.columns = ['KO', 'Reaction']
# print(reaction_ko_df)

### Filter Kegg map file with RR reaction IDs
kri = pd.DataFrame(kegg_reaction_ids, columns=["Reaction"])
# print(kri)
#kri.to_csv("/Users/singh018/Desktop/kegg_RR.tsv", sep="\t")

filtered_df = pd.merge(kri, reaction_ko_df, on="Reaction", how="left")
#ff = filtered_df.dropna()
KO_RR = filtered_df[filtered_df.KO.notnull()]
print(KO_RR)


# filtered_df.to_csv("/Users/singh018/Desktop/filtered.tsv", sep="\t")
#print(ff)

#dup_counts = ff.groupby(ff['Reaction'].tolist(),as_index=False).size()
#dup_counts.to_csv("/Users/singh018/Desktop/final_dup.tsv", sep="\t")
#filtered_df.to_csv("/Users/singh018/Desktop/filtered.tsv", sep="\t")
#ko_RR = temp.loc[temp._merge == 'left_only', ['KO','Reaction']]
#print(ko_RR)


# for each in tables:
#    print(each)
#    columns = cur.execute("SELECT * FROM sqlite_schema WHERE tbl_name = ? and type = 'table'", each)
#    summary[each] = columns
#print(summary)
    # count += 1
# print(count)

# another way
# with sqlite3.connect("mvc.db") as conn:
#    cursor = conn.cursor()

# reading a db in to pandas
# with sqlite3.connect('example.db') as conn:
#    df = pd.read_sql("select * from HOME_PRICES", conn)
#    print(df)

# Writing Pandas DataFrames to DB
# with sqlite3.connect('example.db') as conn:
    # creating a cursor
#    cursor = conn.cursor()
    # reading your own data
#    df = pd.read_csv("your_csv.csv")
    # inserting data
#    rows = [row for name, row in df.iterrows()]
#    cursor.executemany('insert into HOME_PRICES values (?, ?, ?, ?, ?, ?);', rows)

# Creating a db
# with sqlite3.connect('example.db') as conn:
#    cursor = conn.cursor()
#    fields = "id, street, size, number_of_rooms, parking, price"
#    query = f"create TABLE HOME_PRICES ({fields})"
#    cursor.execute(query)
#    houses = [
#        (1, 'first', 110, 4, 1, 100000),
#        (2, 'second', 90, 3, 0, 65000),
#        (3, 'second', 90, 3, 1, 72000)
#    ]
#    cursor.executemany('insert into HOME_PRICES values (?, ?, ?, ?, ?, ?)', houses)
    # cursor.executemany receives an iterable for multiple executions


# A: See https://www.uniprot.org/help/uploadlists, last remark. The IDmapping service favours mappings that are more or less 1:1, not 1:hundreds or even 1:thousands (for performance reasons). You could try to address the issue programmatically, or try queries like *Pfam:PFxxxxx or Pfam:PFyyyyy*. Don’t hesitate to contact us at the helpdesk (https://www.uniprot.org/contact) if you would like to discuss a specific use case, in more detail.
# Example URL: https://www.uniprot.org/uniprot/?query=database:(type:pfam%20pf08563)%20OR%20database:(type:pfam%20pf00870)%20OR%20database:(type:pfam%20pf07710)&format=tab&columns=id,entry%20name,database(Pfam)&sort=score