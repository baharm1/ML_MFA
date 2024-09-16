# Quantifying Cellular Metabolic and Microenvironmental Exchange Fluxes from Gene Expression

We adopted the single-cell flux estimation analysis ([scFEA](https://github.com/changwn/scFEA)) algorithm for serine metabolism and added exchange fluxes of serine within the tumor microenvironment (TME). The intra- and intercellular fluxes are balanced in each cell and in the TME.

## Requirements
The requirements to run `modified_scFEA.py` are similar to <sup>13</sup>C-scMFA.

## Usage

The input arguments of the model are provided to the model via a JSON file (modified_scFEA_config.json), which includes the directories required for the `modified_scFEA.py` to run:

* `data_dir`: Directory of the metabolic model and single-cell data
* `res_dir`: Directory of output files
* `module_gene_file`: Filename of genes associated with metabolic reactions
* `stoichiometry_matrix`: Filename of metabolite balances
* `compound_name_file`: Filename of balanced metabolite names
* `gene_exp`: Filename of a block diagonal gene expression matrix constructed in `prep_input_data_serine_modified_scFEA.R`
* `cellnames_types`: Filename of a dataframe including cell IDs and cell types

The results of modified scFEA are visualized in `visualize_modified_scFEA_output.R`. These plots are shown in Fig. 2F-I for Darmanis _et al._'s scRNA-seq dataset ([GSE84465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465)).