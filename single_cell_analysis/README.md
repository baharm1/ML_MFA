# Single-cell RNA-seq analysis

Here, we analyzed our scRNA-seq data as well as Darmanis _et al._'s scRNA-seq dataset[1] ([GSE84465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465)) to investigate cellular heterogeneity in high grade glioma patients. 
The cell annotations and normalized counts from this step were used in [metabolic interaction analysis](https://github.com/baharm1/ML_MFA/blob/main/single_cell_metabolic_interaction_analysis), [modified scFEA](https://github.com/baharm1/ML_MFA/blob/main/modified_scFEA), and [<sup>13</sup>C-scMFA](https://github.com/baharm1/ML_MFA/blob/main/13C_scMFA).
In addition, single-cell analysis of patient-derived xenografts (GBM12, GBM38, and HF2303) and TRP tumor can be found in this folder. 
The `.gmt` files include cell type and cell state signature genes.

## Requirements
* R 4.2.2
	* dplyr 1.1.2
	* tidyr 1.3.0
	* reshape2 1.4.4
	* ggplot2 3.4.2
	* RColorBrewer 1.1-3
	* scales 1.2.1
	* patchwork 1.1.2
	* Polychrome 1.5.1
	* scCustomize 1.0.0
	* Seurat 4.2.0
	* Matrix 1.5-1
	* DropletUtils 1.18.1
	* DoubletFinder 2.0.3
	* scrabble 1.0.0
	* harmony 0.1.1
	* GSEABase 1.60.0
	* stringr 1.5.0
	* cowplot 1.1.1

## References
[1] Darmanis, S. et al. Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma. _Cell Reports_ (2017). [doi: 10.1016/j.celrep.2017.10.030](https://doi.org:10.1016/j.celrep.2017.10.030)