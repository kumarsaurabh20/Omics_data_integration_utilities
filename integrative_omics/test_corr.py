#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

transcripts = pd.read_csv("./Wisecaver/filtered_transcripts_fpkm.csv")
print("Step 1: Transcripts file is loaded!")
metabolites = pd.read_csv("./Wisecaver/Final_quant_MEJA.csv")
print("Step 1: Metabolites file is loaded!")

trans_t = transcripts.set_index('name').T
metabo_t = metabolites.set_index('ID').T
print(trans_t)
print(metabo_t)
# result = pd.concat([metabo_t, trans_t], axis=1).corr(method='pearson')
#cor = metabo_t.corr(trans_t, method='pearson')


print(cor)

