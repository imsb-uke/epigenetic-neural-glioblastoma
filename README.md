## About

The repository contains code for the following tasks included in [Epigenetic neural glioblastoma enhances synaptic integration and predicts therapeutic vulnerability](https://www.biorxiv.org/content/10.1101/2023.08.04.552017v1) [1].
|Task   | Description             | Source       |
|:--------------------- |:------------------- |:-------------------:|
| DNA methylation deconvolution | Code to process and perform DNA methylation deconvolution. Please cite [Moss et al. (2018)](https://www.nature.com/articles/s41467-018-07466-6) [2] if you use it.| [DNAm_deconv](code/DNAm_deconv)  |
| CNV | Wrapper to perform copy number variation analysis using Conumee package across multiple groups | [CNV](code/CNV) |
| Differential methylation probes | Differential methylation probes and gene set enrichment between the neural groups | [DNAm_DMP](code/DNAm_DMP) |
| Optimal number of clusters | Overcluster the clinical cohort to find if cluster size > 2 is significantly separable with respect to overall survival | [neural_group_over_cluster](code/neural_group_over_cluster) |
| Signature classifying the neural groups | Contains code to stratify neural groups based on DMP probes between low and neural groups. Include trained logistic regression model | [neural_group_classification](code/neural_group_classification) |
| Mutation analysis | Code to generate Oncoplot on mutations | [oncolplot.ipynb](code/Oncoplot.ipynb) |
| RNA groups | Compare correspondence between RNA GBM subgroups and neural subgroups for paired DNAm-RNA TCGA data | [TCGA](code/TCGA) |
| WGCNA | Code to perform WGCNA analysis on paired proteomics data, also includes geneset and cell type enrichment | [WGCNA_proteomics](code/WGCNA_proteomics) |

[1] Drexler, Richard, et al. "Epigenetic neural glioblastoma enhances synaptic integration and predicts therapeutic vulnerability." bioRxiv (2023).

[2] Moss, Joshua, et al. "Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease." Nature communications 9.1 (2018): 5068.
