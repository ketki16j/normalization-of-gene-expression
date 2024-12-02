This repository contains normalization of GTExV8 gene count through different methods. We use methods like DESeq2, inverse quantile normalization to generate normalize gene count.

Gene count normalization is a process that adjusts raw gene expression count values to account for factors other than the RNA expression of interest. This step is essential for making accurate comparisons of gene expression between samples. 

Here are some reasons why gene count normalization is important:
1. Technical factors: Normalization accounts for technical factors that can cause variation in detected gene counts across cells, such as variation in RNA capture rate. 
2. Biological variability: Normalization removes technical variability while keeping biological variability. 
3. Downstream analyses: Gene count normalization is an important step that precedes many downstream analyses. 

Some normalization methods include:
1. Library size normalization: Assumes that all samples have similar total gene expression. 
2. Counts per million (CPM): Normalizes RNA-seq data for sequencing depth but not gene length. 


The dataset for this tutorial is obtained from the GTEx portal: https://www.gtexportal.org/

```Step1: Extract specific tissue gene count:```
We first need to extract gene count specific for a tissue to do normalization. We employe sample INFO/ID to extract specific tissue gene count by utilizing ```GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt``` file and gene count file ```GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct```. My analysis focus on extracting brain cortex tissue gene count for samples that have RNASEQ. 
The script  ```extract_genecount_braincortextissue.R``` extract specific samples for brain cortex that have RNASEQ expression values

```Step2: Normalization by different methods:```
```1. Inverse quantile Normalization:``` We first filter out the genes and samples based on a threshold value to filter out low expressed genes.
```
 >=6 reads in >=20%samples
>=0.1CPM in >=20% samples
```
Then we normalized the gene count by using TMM and apply rank based inverse normal transformation (INT) across samples

```2. DESeq2 normalization: ```
We first need to filter out low expressed genes based on a cutoff:
```
>=6 reads in >=20%samples
>=0.1CPM in >=20% samples
```
Once we filter the low expressed gene , we normalized read count using DESeq2 and apply rank based INT across samples

```Step3: Comparing different normalized expression: ```
The final step is to compare different normalized gene expression by utilizing cor function in R to generate correlation between predicted and normalized actual expression and observe which one gives us a higher accuracy. 
In our case: we clearly see that inverse quantile normalization gives a higer accuracy in comparison to DESeq2 normalization.

![image](https://github.com/user-attachments/assets/8b3bb4d3-10c8-40ce-9258-b684c079ab5a)



