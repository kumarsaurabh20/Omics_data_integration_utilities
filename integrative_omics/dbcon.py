import io
import os
import sqlite3
import pandas as pd

#con = sqlite3.connect("/Users/singh018/Downloads/deseq/Trionate.sqlite")
con = sqlite3.connect("/Users/singh018/Documents/NPDB_merge.sqlite")
cur = con.cursor()
# sql = "SELECT transcript_id, length FROM ORF WHERE transcript_id IS NOT NULL;"
# sql = "SELECT TrinityID, FullAccession FROM BlastDbase WHERE TrinityID IS NOT NULL;"

#sql = "SELECT BlastDbase.TrinityID, BlastDbase.FullAccession, HMMERDbase.HMMERTDomainDescription FROM BlastDbase INNER JOIN HMMERDbase ON HMMERDbase.QueryProtID = ORF.orf_id  INNER JOIN ORF ON ORF.transcript_id = BlastDbase.TrinityID;"
sql = "SELECT ms_name, mz, monoisotopic_mass FROM merged_data;"
cursor = cur.execute(sql)
count = 1
#transcripts = {}
signatures = {}
for each in cursor:
    print(each)
    ms_name, mz, mm = each
    signatures[count] = [ms_name, mz, mm]
    count += 1
con.close()

result = pd.DataFrame(signatures.values(), columns=['ms_name', 'mz', 'mm'])
result.to_csv('mass_signatures.csv', index=False)