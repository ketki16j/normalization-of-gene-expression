#/bin/python

import pandas as pd
import numpy as np
from edgeR import DGEList, calcNormFactors, cpm


# Step 1: Load raw count data
# Replace 'raw_counts.csv' with your file path
raw_counts = pd.read_csv("raw_counts.csv", index_col=0)  # Rows: genes, Columns: samples

# Parameters
min_reads = 6
min_cpm = 0.1
sample_threshold = 0.2  # 20%

# Step 2: Filter genes based on read counts
def filter_genes(raw_counts, min_reads, sample_threshold):
    sufficient_reads = (raw_counts >= min_reads).sum(axis=1) >= sample_threshold * raw_counts.shape[1]
    return raw_counts[sufficient_reads]

filtered_counts = filter_genes(raw_counts, min_reads, sample_threshold)

# Step 3: Filter genes based on CPM
dge = DGEList(counts=filtered_counts)
dge = calcNormFactors(dge)  # TMM normalization
cpm_data = cpm(dge)  # Convert to CPM
sufficient_cpm = (cpm_data >= min_cpm).sum(axis=1) >= sample_threshold * cpm_data.shape[1]
filtered_cpm = cpm_data[sufficient_cpm]

# Step 4: Apply rank-based inverse normal transformation (INT)
def rank_int_transform(df):
    def int_transform(x):
        ranks = rankdata(x)
        return norm.ppf((ranks - 0.5) / len(ranks))
    return df.apply(int_transform, axis=1)

normalized_counts = rank_int_transform(filtered_cpm)
