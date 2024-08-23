
# Simulation of fluxes and MIDs

Due to the limitation of <sup>13</sup>C-metabolomics data of HGG patients and the uncertainties around quantifying metabolic fluxes in patients, training a digital twin framework on actual patient data is not feasible. Here, we employ flux balance analysis (FBA) and mass isotopomer/isotopologue balance analysis at isotopic unsteady state (INST-MFA) to generate thousands of synthetic patient data instances. We verify that the simulated MIDs are representative of actual patient MIDs. Then, the simulated MIDs and fluxes are used as input features and targets, respectively, in the metabolic CNN model for training, validating, and testing the model.

Here, we show MID and flux simulations for two pathways: serine metabolism and purine metabolism. Furthermore, to evaluate the simulated MIDs and fluxes, the simulated MIDs are compared to the experimental MIDs.

## Requirements

MATLAB R2021b with the Parallel Processing Toolkit was used to simulate MIDs and fluxes. To visualize the simulated MIDs and fluxes, R version 4.2.2 and the following R packages were used:

* ggplot2 3.4.2
* stringr 1.5.0
* gghalves 0.1.4
* tidyr 1.3.0
* dplyr 1.1.2
* see 0.8.0
* Seurat 4.2.0

## User-defined Inputs for MID and Flux Simulation

* `isotopomer_model`: An Excel file containing the list of reactions, types of reactions, carbon transitions, number of carbons in metabolites, input and sink metabolites, and unlabeled metabolites.
* `isotopologue_model`: An Excel file containing the list of reactions, types of reactions, and input metabolites.
* `param_bound_file`: An Excel file containing the list of reactions, and the lower and upper bounds of fluxes and balanced metabolite pools.
* `F`: Scaling factor for exchange flux based on Wiechert's bidirectional modeling.
* `n_sim`: Number of flux sets or simulations.
* `tspan`: Time period for MID simulation with a specified time step.
* Lower and upper bounds for input metabolite MIDs.
* Source fluxes of a target metabolite.
* Target flux. 

### Isotopomer Model

The isotopomer model is an Excel file with six sheets, described as follows. Isotopomer balance is used for a metabolic model with cleavage reactions.

#### `rxns`

This sheet includes a list of model reactions, reaction types, and carbon transitions as shown in the table below. Reactions are classified as irreversible (I), reversible (R), sink or output (E), and rapid conversion (Q).

| Reaction Type | Reaction   | Carbon Transition |
| ------------- | ---------- | ----------------- |
| Q or I or R   | X == Y + Z | abc,ab,c          |
| E             | X == 0     |                   |

#### `metab_size`

This sheet includes all metabolites used in the model and the number of carbon atoms in metabolites.

| Metabolite | Number of Carbons |
| ---------- | ----------------- |
| X          | 3                 |
| Y          | 2                 |
| Z          | 1                 |

#### `input_metab`

A column that includes the names of input metabolites. Input metabolites are outside of the model boundaries and are not stoichiometrically or isotopically balanced.

#### `unlabeled`

A column that includes the names of unlabeled metabolites. Some metabolites have zero or negligible isotopic enrichment and are considered unlabeled.

#### `EQ`

A column that includes the names of substrates with rapid conversion. In the rapid conversion of substrates to products, the reaction rate is high, and it is assumed that the enrichment of substrates and products is the same.

#### `remove`

A column that includes metabolites that have large pools and participate in many reactions, such as CO2. These metabolites are not stoichiometrically or isotopically balanced.

### Isotopologue Model

The isotopologue model is an Excel file with two sheets, described as follows. Isotopologue balance can be used when there is no cleavage reaction in the metabolic model.

#### `rxns`

This sheet includes a list of model reactions and reaction types. Reactions are classified as irreversible (I) or sink (E). These reactions can be condensation or rearrangement reactions. A reversible reaction can be included in this model if the forward or backward reaction is not a cleavage reaction. If so, the forward and backward reactions are added separately to the table below.

| Reaction Type | Reaction     |
| ------------- | ------------ |
| I             | X + Y == Z   |
| I (forward)   | X == W       |
| I (backward)  | W == X       |
| E             | X == 0       |

#### `input_metab`

A column that includes the names of input metabolites. Input metabolites are outside of the model boundaries and are not stoichiometrically or isotopically balanced.

## General Usage

By providing the input arguments to our simulation model, MIDs and fluxes are generated. In general, one can curate their own metabolic model based on the metabolic pathways presented in the literature. Following the formatting provided in the previous section, the model needs to be in Excel format with the mentioned sheet names. If there is not any cleavage reaction in the model, it is more efficient to use isotopologue balances instead of isotopomer balances. These two models are separated in `isotopomer_model` and `isotopologue_model` files. For isotopomer balances, we used atom mapping (AMM) and isotopomer mapping (IMM) methods to generate a system of ordinary differential equations (ODEs) for isotopomer gradients. The isotopomer balances are automatically generated by our function `ODE_fn_gen` if the carbon atom transition is defined in the isotopomer model. One example of an isotopomer model is described below in the serine model. For the isotopologue model, one should write these balance equations in the `simulate` function. An example of an isotopologue model is provided in the purine model section.

## Replication of Serine Simulation in the Manuscript

We applied isotopomer balances to the serine metabolic model since it includes the cleavage of serine to glycine and 5, 10-methylenetetrahydrofolate (MTHF), and the degradation of glycine through the glycine cleavage system. Since carbon atom transition is important in the cleavage reactions, isotopomer balances are required. By providing a list of reactions and carbon transitions from substrates to products in `model_serine.xlsx`, `import_stoich_AM` creates stoichiometry and atom mapping matrices that are further used in `ODE_fn_gen` to generate isotopomer balance equations at the isotopic unsteady state (`serine_ODE`). Then, the lower and upper bounds of metabolic fluxes and metabolite concentrations are imported from `parameter_bounds_serine_model.xlsx`. Fluxes and concentrations are generated randomly within the bounds. To generate target fluxes with a uniform distribution that creates a balanced dataset for the metabolic CNN to train, the ratio of a target flux (i.e., glucose-derived serine synthesis, plasma serine uptake, or tumor microenvironment (TME)-derived uptake) to all source fluxes of serine is set to a linearly spaced vector between 0.01 and 0.99. Furthermore, fluxes are balanced given stoichiometric coefficients within the flux bounds through a constrained optimization problem. 

To simulate serine model MIDs, the lower and upper bounds of fractional enrichment of input metabolite MIDs are defined by the user based on the ranges observed in patients. Input metabolite MIDs are generated randomly within the bounds for various simulations and passed to the `simulate` function along with balanced fluxes and metabolite concentrations. In the `simulate` function, a system of ODEs of isotopomer balances (`serine_ODE`) is solved, which estimates the gradient of isotopomers of glycine and MTHF over time. Finally, isotopomer distributions map into isotopologue distributions and isotopologues of balanced metabolites are saved at different time points. Please refer to our manuscript for more details on flux and MID simulation.

#### Input Arguments

* `isotopomer_model`: `model_serine.xlsx`
* `param_bound_file`: `parameter_bounds_serine_model.xlsx`
* `F`: 200 pmol/mg/h
* `n_sim`: 50,000
* `tspan`: 0 to 4 h with a step size of 0.1 h
* Lower and upper bounds for input metabolite MIDs: Plasma serine (SERp), cortex phosphoglycerate (PGc), and glioma phosphoglycerate (PGg) are the input metabolites of the serine model, and the bounds of fractional enrichment of their isotopologues are set as follows according to our patient MIDs:
| Metabolite | M + 1 Bounds | M + 2 Bounds | M + 3 Bounds |
| ---------- | ------------ | ------------ | ------------ |
| PGg        | [0, 0.06]    | [0, 0.035]   | [0.02, 0.17] |
| PGc        | [0, 0.03]    | [0, 0.015]   | [0.02, 0.25] |
| SERp       | [0.2, 0.7]   | [0, 0.02]    | [0, 0.01]    |
* Source fluxes of a target metabolite: For glioma serine, three labeled sources are glucose-derived serine synthesis (`PGg == SERg`), plasma serine uptake (`SERp == SERg`), and TME-derived serine uptake (`SERc == SERg`). For cortex serine, two labeled sources are glucose-derived phosphoglycerate (`PGc == SERc`) and plasma serine uptake (`SERp == SERc`).
* Target flux: Any of the labeled sources described above. A simulation can be performed on one target flux at a time.

## Replication of Purine Simulation in the Manuscript

We used part of our serine metabolic model to simulate glycine and MTHF MIDs which are further used in inosine monophosphate (IMP) de novo synthesis. As described above, serine metabolic model includes cleavage reactions and we need to apply isotopomer balances. However, isotopologue balances can be applied to purine reactions. Fluxes and concentrations are generated randomly within the bounds defined in `parameter_bounds_gmp_denovo_glioma`. To generate target fluxes with a uniform distribution that creates a balanced dataset for the metabolic CNN to train, the ratio of a target flux (i.e., de novo or salvage guanosine monophosphate (GMP) synthesis) to all source fluxes of GMP is set to a linearly spaced vector between 0.01 and 0.99. Furthermore, fluxes are balanced given stoichiometric coefficients within the flux bounds through a constrained optimization problem. 

To simulate purine model MIDs, the lower and upper bounds of fractional enrichment of input metabolite MIDs are defined by the user based on the ranges observed in patients. Input metabolite MIDs are generated randomly within the bounds for various simulations and passed to the `simulate` function along with balanced fluxes and metabolite concentrations. In the `simulate` function, a system of ODEs of isotopologue balances is solved, which estimates the gradient of isotopomers and isotopologues of glycine, MTHF, IMP, GMP, guanosine diphosphate (GDP), guanosine, inosine, and adenosine monophosphate (AMP) over time. Consequently, isotopologues of balanced metabolites are saved at different time points. Please refer to our manuscript for more details on flux and MID simulation.

#### Input Arguments

* `isotopomer_model`: `model_serine.xlsx`
* `isotopologue_model`: `model_purine.xlsx`
* `param_bound_file`: `parameter_bounds_purine_model.xlsx`
* `F`: 200 pmol/mg/h
* `n_sim`: 50,000
* `tspan`: 0 to 4 h with a step size of 0.1 h
* Lower and upper bounds for input metabolite MIDs: Serine (SER), and ribose 5-phosphate (R5P) are the input metabolites of the purine model, and the bounds of fractional enrichment of their isotopologues are set as follows according to our patient MIDs:
| Metabolite | M + 1 Bounds | M + 2 Bounds | M + 3 Bounds | M + 4 Bounds | M + 5 Bounds |
| ---------- | ------------ | ------------ | ------------ | ------------ | ------------ |
| SER        | [0.04, 0.2]  | [0, 0.08]    | [0, 0.12]    |              |              |
| R5P        | [0, 0.06]    | [0, 0.08]    | [0, 0.12]    | [0, 0.03]    | [0, 0.08]    |
* Source fluxes of a target metabolite: For GMP, two labeled sources are de novo GMP synthesis (`IMP == GMP`) and salvage GMP pathway (`GUA + R5P == GMP`). 
* Target flux: Any of the labeled sources described above. A simulation can be performed on one target flux at a time.

## Output Files of Flux and MID simulations

* Simulation parameters: A CSV file containing balanced fluxes, concentrations of balanced metabolites, and MIDs of input metabolites.
* Simulated MIDs: A CSV file containing simulation indices, time points, and MIDs of balanced metabolites.
These files can be found in `../metabolic_CNN/data/sim_data` folder for serine and purine simulations.

## Evaluation of simulated data

To investigate whether simulated MIDs resemble patient MIDs:
* Ranges of simulated MIDs overlap with patient MIDs.
* T-SNE analysis of simualted MIDs and patient MIDs projects them into a same space.
These analyses can be found in `evaluate_simulated_data.R`.
