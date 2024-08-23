# Quantifying Cellular Metabolic, Circulating, and Microenvironmental Exchange Fluxes from Gene Expression and Bulk <sup>13</sup>C-Metabolomics

By integrating scRNA-seq data with mass isotopologue distribution (MID) of cortex, tumor, and circulating metabolites, <sup>13</sup>C-scMFA quantifies metabolic fluxes, circulating, and microenvironmental exchange fluxes for each cell. Based on the location of cells, whether they belong to normal or tumor tissues, bulk MIDs are assigned to each cell. MIDs of microenvironment-derived metabolites are assumed to be similar to those of the tissues that secrete the metabolites. These secretion or uptake fluxes have been defined in modified scFEA analysis. Circulating metabolite MIDs help distinguish between the two exchange fluxes: microenvironment-derived metabolite exchange between cell types and circulating metabolite uptake. The intra- and intercellular fluxes are balanced in each cell and in the tumor microenvironment, similar to modified scFEA. Furthermore, fluxes are adjusted such that the accumulation of isotopologues in cells and tissues is minimized. Here, two examples of using <sup>13</sup>C-scMFA are shown for serine and purine metabolism. Please refer to our manuscript for more information about the methods.

## Requirements

* MATLAB R2021b for estimation of 5,10-methylenetetrahydrofolate (MTHF) MIDs with the following package:
	* Artelys Knitro software 12.4
* R 4.2.2 for visualization and data preparation with the following packages:
	* Seurat 4.2.0
	* dplyr 1.1.2
	* ggplot2 3.4.2
	* RColorBrewer 1.1-3
	* scCustomize 1.0.0

## Usage

1. Curation of metabolic models: We first curated metabolic models of serine and purine pathways. These model are located at `serine_model` and `purine_model` folders. To establish these models, the following files are needed:
	* `cmMAT` defines whether a stoichiometry-balanced metabolite is consumed or produced in a reaction.
	* `cName` includes names of stoichiometry-balanced metabolites or rownames of cmMat.
	* `mmMAT` defines whether an isotopologue of a balanced metabolite is consumed or produced in a reaction.
	* `mName` includes names of isotopologues of balanced metabolites.
	* `module_gene` includes enzymes that catalyzes reactions in the model.
1. Input data prep: The scRNA-seq and <sup>13</sup>C-isotope tracing data were prepared for <sup>13</sup>C-scMFA. These data can be found in `patient_input_serine` and `patient_input_purine`.
	* We collected percent enrichment of isotopologues for input and balanced metabolites in `patient_MID` folder.
	* Preparation of input scRNA-seq data to align with <sup>13</sup>C-scMFA is provided in `prep_input_data_scMFA.R`. By running this script, normalized block diagonal scRNA-seq data for each patient and cell types are saved in `patient_scRNA` and `patient_celltypes`.
1. Running <sup>13</sup>C-scMFA: The input arguments of the model are provided to the model via a JSON file (serine_13C_scMFA_config.json for serine, purine_13C_scMFA_config.json for the purine model), which includes the directories required for the `13C_scMFA.py` to run:
	* `model_name`: Name of the metabolic model (serine or purine)
	* `model_dir`: Directory of the metabolic model
	* `patient_dir`: Directory of patient data, including patient MIDs (percent enrichment), patient scRNA-seq data in the form of a block diagonal gene expression matrix constructed in `prep_input_data_serine_scMFA.R`, and cell types assigned to cell IDs in patient scRNA-seq data. These subdirectories are as follows: `patient_MID`, `patient_scRNA`, and `patient_celltypes`.
	* `output_dir`: Directory of <sup>13</sup>C-scMFA output files
	* `module_gene_file`: Filename of genes associated with metabolic reactions in `model_dir`
	* `stoichiometry_matrix`: Filename of metabolite balances in `model_dir`
	* `isotopologue_matrix`: Filename of isotopologue balances in `model_dir`
	* `compound_name_file`: Filename of balanced metabolite names in `model_dir`
	* `isotopologue_name_file`: Filename of balanced isotopologue names in `model_dir`
	* `n_mids`: Number of balanced isotopologues for a balanced metabolite
	* `n_balanced_metabs`: Number of balanced metabolites in a cell
1. Visualizing <sup>13</sup>C-scMFA estimated fluxes: `visualize_scMFA_output.R` generates Fig. 6B-F, Extended Fig. 12, 13, Fig. 7B, C, Extended Fig. 15 C-J.

#### Metabolite balances

Metabolite balances provided as a CSV file in `stoichiometry_matrix` are stoichiometry coefficients of balanced metabolites in all reactions. This matrix has a shape of the number of balanced metabolites by the number of fluxes. If a balanced metabolite is consumed in a reaction, the sign of its stoichiometric coefficient is negative, and if it is produced, the sign is positive.

#### Isotopologue balances

Isotopologue balances are provided as a CSV file in `isotopologue_matrix` with a shape of the number of MIDs of balanced metabolites by the number of fluxes. This matrix includes coefficients of fluxes in isotopologue balance equations. If a metabolite is consumed in a reaction, the isotopologue coefficient is prefixed with a `neg_` string to denote the negative sign. In purine model, labeled substrates of de novo IMP synthesis are one ribose 5-phosphate, one glycine, and two MTHF. We estimated MTHF enrichment from serine and glycine enrichment data since MTHF enrichment was not available for our patients. We calculated multiplication of isotoplogues of substrates of de novo (denoted by DN in patient MID files) IMP synthesis in `patient_input_purine/MTHF_estimation/denovo_IMP_MID.xlsx`. Please refer to our manuscript for more details.

