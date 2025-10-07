# Single Cell Analysis of High Grade Glioma patients 

Here, we analyzed our scRNA-seq data as well as Darmanis _et al._'s scRNA-seq dataset[1] ([GSE84465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465)) to investigate cellular heterogeneity in high grade glioma patients. 
Furthermore, Darmanis _et al._'s dataset used in [metabolic interaction analysis](https://github.com/baharm1/ML_MFA/blob/main/single_cell_metabolic_interaction_analysis) and [modified scFEA](https://github.com/baharm1/ML_MFA/blob/main/modified_scFEA) to quantify serine exchange rate between cell types. 
These datasets were used in [<sup>13</sup>C-scMFA](https://github.com/baharm1/ML_MFA/blob/main/13C_scMFA) to estimate intra- and intercellular serine fluxes in astrocytes, neurons, and neoplastic cells and purine single-cell fluxes in myeloid and neoplastic cells.
In addition, single-cell analysis of patient-derived xengrafts (GBM12, GBM38, and HF2303) and TRP tumor can be found in this folder.

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

[1] Darmanis, S. et al. Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma. _Cell Reports_ (2017). [doi: 10.1016/j.celrep.2017.10.030](https://doi.org:10.1016/j.celrep.2017.10.030)