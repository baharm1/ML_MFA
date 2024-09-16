# Metabolic flux analysis informed by CNN predicted relative fluxes of serine model

A metabolic flux analysis (MFA) was applied on _in vivo_ data to evaluate CNN by back calculation of MIDs from CNN predicted fluxes and comparing them with experimental MIDs. In addition, we are able to compare serine fluxes between two compartments (_i.e._ glioma and cortex). Since CNN predicts fluxes relatively in each compartment, comparing a similar flux but from two compartments is not possible. In CNN and in our previous MFA model [1], fluxes of serine sources were added to 1 in each compartment, hence comparing fluxes in the same compartment was achievable. It is good to note that an MFA alone is highly underdetermined, and many solutions are possible. In addition, MFA assumes isotopic steady state opposed to CNN framework.
The MFA informed by CNN model comprises two compartments cortex and glioma. Reactions used in MFA model are defined in `two_comp_model.xlsx`. Both compartments have their own phosphoglycerate pool to synthesize serine _de novo_. A circulating serine pool is also available for both compartments. Based on our `modified_scFEA` analysis which shows there is a net flux from cortex to glioma, we also added a TME-derived serine pool for glioma compartment. An unlabeled serine source was also added to both compartments to include serine sources from protein breakdown and autophagy. Since serine model includes cleavage reactions, MFA balance equations are written on isotopomers which results in introducing more variables than equations. To reduce the model uncertainty, we added 4 linear constraints:
* Fluxes of serine sources in glioma are set to CNN predicted relative fluxes (sum of serine sources equals 1). This adds 3 linear equality constraints for _de novo_ serine synthesis, plasma serine uptake, and TME-derived serine uptake.
* Ratio of de novo serine synthesis to total serine synthesis in cortex is set to CNN prediction. This adds 1 linear equality constraint.

Model parameters consist of fluxes and isotopomers of plasma serine and phosphoglycerate, serine, glycine, and 5,10-methylenetetrahydrofolate in glioma and cortex. For the metabolites inside the compartments, we balanced fluxes by stoichiometry coefficients assuming that fluxes are at steady state. Isotopomers of metabolites inside the compartments are also balanced assuming that they are at steady state. The model parameters are initialized randomly 10 times within the range of 0 and 1. We set the initial values of unlabeled serine uptake fluxes to 0 to prevent overestimation of unlabeled pools. Parameters are optimized to minimize the differences between experimentally measured MIDs and MIDs estimated by MFA.

## Requirements
* MATLAB R2021b with the following package:
	* Artelys Knitro software 12.4

## Usage 
The following input files are required to run `MFA_with_ML_fluxes.m`:
* Metabolic model: an isotopomer model is defined for serine reactions in `two_comp_model.xlsx`. The description of sheets can be found in `data_simulation`.
* Patient MID data: MIDs of serine model metabolites can be found in `patient_input_data`. Columns of excel files represent metabolite and site of sample, isotopologue, mean, and standard deviation of technical replicates.
* CNN prediction of relative fluxes of serine sources: average of relative fluxes predicted for each patient is given to MFA as the linear constraint.

The results of `MFA_with_ML_fluxes.m` are saved in `output_files`.
Given these output files to `sensitivity_analysis.m`, confidence intervals of fluxes can be determined.

## References
[1] Scott, A. J. et al. Rewiring of cortical glucose metabolism fuels human brain cancer growth. medRxiv (2023).