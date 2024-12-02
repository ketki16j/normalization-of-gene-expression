GTEx_Analysis_v8_Annotations_SampleAttributesDS <- read_delim("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",delim = "\t", escape_double = FALSE,trim_ws = TRUE)
sample_attributes <- select(GTEx_Analysis_v8_Annotations_SampleAttributesDS,SAMPID,SMTS,SMTSD,SMAFRZE)
#head(sample_attributes)

sample_attributes_braindata <- sample_attributes %>% filter(sample_attributes$SMTS == "Brain" & sample_attributes$SMAFRZE == "RNASEQ")
#head(sample_attributes_braindata)
#dim(sample_attributes_braindata)
sample_attributes_braindata <- data.frame(sample_attributes_braindata)
head(sample_attributes_braindata)[1:5,1:4]


GTEx_Analysis_gene_reads <- read_table2("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
#head(GTEx_Analysis_gene_reads)
#str(GTEx_Analysis_gene_reads)
GTEx_Analysis_gene_reads <- data.frame(GTEx_Analysis_gene_reads)
gene_count <- data.frame(t(GTEx_Analysis_gene_reads[,-c(1:2)]))
colnames(gene_count) <- GTEx_Analysis_gene_reads$Name
head(gene_count)[1:5,1:5]
library(tidyverse)
gene_count <- rownames_to_column(gene_count,var = "SAMPID")
gene_count <- mutate(gene_count,SAMPID = gsub(pattern = ".", replacement = "-", x = SAMPID, fixed = TRUE))
brain_tissueid <- select(sample_attributes_braindata,SAMPID)

new_dataset_gene_Count_braintissue <- gene_count %>% right_join(brain_tissueid, by="SAMPID")
write.table(new_dataset_gene_Count_braintissue,file="~/Desktop/genecount_braintissue_Exprdata2.txt", sep= "\t", quote = F, row.names = F)


brain=read.delim("genecount_braintissue_Exprdata2.txt", head = T, row.names = 1)
dim(brain)





