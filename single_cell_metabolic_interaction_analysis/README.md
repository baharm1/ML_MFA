# Single Cell Metabolic Interaction Analysis in GBM Patients
To study metabolic interactions between cell types in the glioblastoma (GBM) microenvironment, we selected Darmanis _et al._'s scRNA-seq dataset [1] ([GSE84465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465)) since it includes both tumor and periphery samples. We applied MEBOCOST [2] on Darmanis _et al._'s dataset to investigate the cell-cell metabolic inference between different cell types in the GBM microenvironment. In addition, we defined a metabolic accumulation score (MAS) to identify cell types that produce/consume a metabolite more than other cell types in the microenvironment. Based on MEBOCOST and MAS analyses, we identified astrocytes/neurons as potential cell types that secrete/uptake serine into/from the tumor microenvironment.

## Usage
First step is to install MEBOCOST from [its GitHub repository](https://github.com/zhengrongbin/MEBOCOST). Navigate to this folder to install MEBOCOST in the current directory, then run `MEBOCOST_communication_score.ipynb` which utilizes MEBOCOST for communication prediction. 
This jupyter notebook requires three input files, as follows:
* `counts_darmanis.h5`: log-normalized reads with scale factor of 10,000 which is equivalent to data slot of a Seurat object
* `metadata_damanis.csv`: a dataframe with number of cells rows and attributes as columns including cell types
* `mebocost.conf`: a configuration file which refers to `my_config_files` for metabolic model 
The normalized counts and metadata of Darmanis dataset were saved previously by running `../single_cell_analysis/scRNA_analysis_Darmanis.R`. The plots are saved in `darmanis_MEBOCOST_output` folder.

Calculation of MAS can be found in `calc_metabolite_accumulation_score_serine.R`. This code requires a preprocessed Seurat object, herein, Darmanis scRNA-seq data.

## References
[1] Darmanis, S. et al. Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma. Cell Reports 21, 1399-1410 (2017). [doi: 10.1016/j.celrep.2017.10.030](https://doi.org:10.1016/j.celrep.2017.10.030)

[2] Zheng, R. et al. MEBOCOST:Metabolite-mediated Cell Communication Modeling by Single Cell Transcreptome. bioRxiv (2022). [doi: https://doi.org/10.1101/2022.05.30.494067](https://doi.org/10.1101/2022.05.30.494067)