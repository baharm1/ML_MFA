# Quantifying Cellular Metabolic and Microenvironmental Exchange Fluxes from Gene Expression

We modified the single-cell flux estimation analysis algorithm[1] ([scFEA](https://github.com/changwn/scFEA)) to quantify cellular serine metabolic reaction and serine exchange fluxes within the tumor microenvironment (TME) by minimizing the accumulation of serine in cells and TME.

## Requirements
The requirements to run `modified_scFEA.py` are similar to <sup>13</sup>C-scMFA.

## Usage

The input arguments of the model are provided to the model via a JSON file (modified_scFEA_config.json), which includes the directories required for the `modified_scFEA.py` to run:

* `data_dir`: Directory of the metabolic model and single-cell RNA-seq data
* `res_dir`: Directory of output files
* `module_gene_file`: Filename of genes associated with metabolic reactions
* `stoichiometry_matrix`: Filename of metabolite balances
* `compound_name_file`: Filename of balanced metabolite names
* `gene_exp`: Filename of a block diagonal gene expression matrix constructed in `prep_input_data_serine_modified_scFEA.R`. The gene expression of Darmanis dataset can be found in `data/three_comp_scRNA.csv`.
* `cellnames_types`: Filename of a dataframe including cell IDs and cell types. Cell types of Darmanis dataset can be found in `data/cellnames_types`

The results of modified scFEA (`output`) are visualized in `visualize_modified_scFEA_output.R`. These plots are shown in **Figure 2D-G** for Darmanis _et al._'s scRNA-seq dataset. Requirements of each code can be found in the scripts.

## References
[1] Alghamdi N. _et al._ A graph neural network model to estimate cell-wise metabolic flux using single-cell RNA-seq data, _Genome Research_ (2021)
