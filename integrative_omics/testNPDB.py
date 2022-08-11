import os
import sys
import argparse
import sqlite3
import pandas as pd
from multiprocessing import Pool

import gizmos

conn = sqlite3.connect("NPDB.sqlite")
c = conn.cursor()
query = c.execute("SELECT structure.structure_id,structure.monoisotopic_mass, structure_has_data_source.source_id, structure_has_data_source.source_name, structure.inchi,structure.inchi_key2,structure.smile FROM structure left join structure_has_data_source on structure_has_data_source.structure_id = structure.structure_id")
cols = [column[0] for column in query.description]
print(cols)
results= pd.DataFrame.from_records(data = query.fetchall(), columns = cols)
print(results)
c.close()