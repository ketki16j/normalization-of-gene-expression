# Load necessary libraries
library(DESeq2)
library(edgeR)

# Assume `count_data` is your raw count matrix (genes as rows, samples as columns)
# `metadata` is the metadata for the samples

# Step 1: Filter genes based on read count threshold
filter_genes <- function(count_data, threshold = 6, min_samples = 0.2) {
  sample_count <- ncol(count_data)
  min_samples <- ceiling(min_samples * sample_count)
  filtered <- count_data[rowSums(count_data >= threshold) >= min_samples, ]
  return(filtered)
}

# Step 2: Filter genes based on CPM threshold
filter_genes_cpm <- function(count_data, cpm_threshold = 0.1, min_samples = 0.2) {
  sample_count <- ncol(count_data)
  min_samples <- ceiling(min_samples * sample_count)
  cpm_values <- cpm(count_data)
  filtered <- count_data[rowSums(cpm_values >= cpm_threshold) >= min_samples, ]
  return(filtered)
}

# Apply both filters
filtered_counts <- filter_genes(count_data)
filtered_counts <- filter_genes_cpm(filtered_counts)

# Step 3: Normalize using DESeq2
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = metadata,
                              design = ~ 1)  # Assuming no specific design for normalization
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Step 4: Apply rank-based inverse normal transformation (INT)
rank_int <- function(x) {
  n <- length(x)
  ranks <- rank(x, ties.method = "average")
  qnorm((ranks - 0.5) / n)
}

int_normalized_counts <- apply(normalized_counts, 1, rank_int)
int_normalized_counts <- t(int_normalized_counts)  # Transpose back to original format

# Output the processed data
write.csv(int_normalized_counts, "int_normalized_counts.csv", row.names = TRUE)

# Verify results
print("Filtering and normalization completed successfully.")
