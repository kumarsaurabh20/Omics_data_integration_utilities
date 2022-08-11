import glob
import os
import time
import pandas as pd


def sort_columns(file):
	# cols = ["transcript_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
	temp = pd.read_csv(file)
	temp = temp.sort_values('padj')
	temp = temp[temp['padj'] <= 0.05]
	return temp

def main():
	file_subset = "/Users/singh018/Documents/Meantools_v1/Wisecaver/Wisecaver/metabolome_MEJA.csv"
	file_all = "/Users/singh018/Documents/Meantools_v1/Wisecaver/Wisecaver/mass_signatures.sorted.csv"

	file_sub_df = pd.read_csv(file_subset, names=["ms_name", "mz"])
	print("STEP 1: Created a dataframe for blast db!")
	print(file_sub_df)

	file_all_df = pd.read_csv(file_all, names=["ms_name", "mz", "mm"])
	print("STEP 2:  Created a dataframe for uniprot-pfam maps!")
	print(file_all_df)

	result = pd.merge(left=file_sub_df, right=file_all_df, on="ms_name", how="inner")
	print("STEP 3: Merged transcripts-uniprot and uniprot_pfam dataframe!")
	print(result)
	result.to_csv("merged_metabolome_MeJA.csv", index=False)


if __name__ == "__main__":
	start_time = time.time()
	main()
	end_time = time.time()
	print(end_time - start_time)
