# fix the import file path to be the feature counts output, the directory prefix to be the longer names, and the output path
## specifically, lines 7, 8 and 35 in this file
# must have the DESeq2 library installed

library("DESeq2")

all_counts <- read.delim("/home/scratch/dir/GVDS_all_bam_counts_2021v1.txt", header = TRUE, skip=1, sep="\t")
dir_prefix<-"X.scratch.user.GEUVADIS.bam_input."

# removes the columns that are not the geneid and sample headers
all_counts<-all_counts[,c(1, 7:ncol(all_counts))]

# auto sample name does like a full directory; ex: â€œX.home.user.scratch.TWAS_test_bam.HG02215.1.M_111124_4.bamâ€
# fix the sample names to be shorter
names(all_counts)<-gsub(dir_prefix, "", colnames(all_counts))
names(all_counts)<-gsub(".bam", "", colnames(all_counts))

# changes the row names to gene names and then updates the columns to only be the samples
# seems like the dataframe will work best with this set up
geneNames<-all_counts[,1]
rownames(all_counts)<-geneNames
all_counts<-all_counts[,2:ncol(all_counts)]

# creates a reference table to run the DE seq with
sample_info<-DataFrame(condition=names(all_counts), row.names=names(all_counts))

# runs the DESeq2
ds<-DESeqDataSetFromMatrix(countData=all_counts, colData=sample_info, design= ~condition)
keep_genes<-rowSums(counts(ds))>0
# from 63677 to 56957
ds<-ds[keep_genes,]
ds<-estimateSizeFactors(ds)
normalized_counts<-counts(ds, normalized=TRUE)

write.table(normalized_counts, file = "GVDS_normalized_counts.txt", sep="\t", quote=FALSE, col.names=NA)

