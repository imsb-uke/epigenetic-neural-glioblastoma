## About

The repository contains code for the following tasks included in [A prognostic neural epigenetic signature in high-grade glioma](https://www.nature.com/articles/s41591-024-02969-w) [1]:
|Task   | Description             | Source       |
|:--------------------- |:------------------- |:-------------------:|
| An example run of the signature classifying the neural groups | Contains code to run the fitted logistic regression model on an example IDAT | [run neural classification](code/example_neural_classification) |
| DNA methylation deconvolution | Code to process and perform DNA methylation deconvolution. Please cite [Moss et al. (2018)](https://www.nature.com/articles/s41467-018-07466-6) [2] if you use it.| [DNAm_deconv](code/DNAm_deconv)  |
| CNV | Wrapper to perform copy number variation analysis using Conumee package across multiple groups | [CNV](code/CNV) |
| Differential methylation probes | Differential methylation probes and gene set enrichment between the neural groups | [DNAm_DMP](code/DNAm_DMP) |
| Optimal number of clusters | Overcluster the clinical cohort to find if cluster size > 2 is significantly separable with respect to overall survival | [neural_group_over_cluster](code/neural_group_over_cluster) |
| Signature classifying the neural groups | Contains code to stratify neural groups based on DMP probes between low and neural groups. Include trained logistic regression model | [neural_group_classification](code/neural_group_classification) |
| Mutation analysis | Code to generate oncoplot on mutations | [Oncolplot.ipynb](code/Oncoplot.ipynb) |
| RNA groups | Compare correspondence between RNA GBM subgroups and neural subgroups for paired DNAm-RNA TCGA data | [TCGA](code/TCGA) |
| WGCNA | Code to perform WGCNA analysis on paired proteomics data, also includes geneset and cell type enrichment | [WGCNA_proteomics](code/WGCNA_proteomics) |

## References
[1] Drexler, R., Khatri, R., Sauvigny, T. et al. A prognostic neural epigenetic signature in high-grade glioma. Nat Med (2024). https://doi.org/10.1038/s41591-024-02969-w

[2] Moss, J., Magenheim, J., Neiman, D. et al. Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease. Nat Commun 9, 5068 (2018). https://doi.org/10.1038/s41467-018-07466-6
