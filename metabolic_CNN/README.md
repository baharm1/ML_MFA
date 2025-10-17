# Metabolic CNN: Estimation of relative anabolic fluxes in bulk tissues

## Requirements

* Python 3.8 for running metabolic CNN: We recommend installing [Anaconda](https://www.anaconda.com/download) and creating a new conda environment to manage the versioning of dependencies:
```
conda create -n mCNN python=3.8.15
conda activate mCNN
```
To run Python codes, the following packages have to be installed:
```
pip install -r requirements.txt
```
To install tensorflow, please visit [tensorflow documentation](https://www.tensorflow.org/install/pip). 
We installed tensorflow on a GPU-enabled workstation (NVIDIA Quadro P2200, OS:Windows).
```
conda install -c conda-forge cudatoolkit=11.2 cudnn=8.1.0
pip install tensorflow == 2.10.1
```
To verify tensorflow installation:
```
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
```

* R 4.2.2 for visualization and data preparation with the following packages:
	* Seurat 4.2.0
	* dplyr 1.1.2
	* tidyr 1.3.0
	* ggplot2 3.4.2
	* ggpubr 0.4.0
	* ggbeeswarm 0.7.2
	* stringr 1.5.0

## Usage
1. Create config file: The input arguments are provided to the model via a JSON file, which includes the following dictionary required for `run_metabolic_CNN.py` to run:
	* `sim_data_dir`: Directory of simulated data that contains sim_data_filename and sim_param_filename
	* `sim_data_filename`: Filename of simulated data that contains MIDs of balanced metabolites
	* `sim_param_filename`: Filename of simulated parameters that contains MIDs of input metabolites and fluxes
	* `cols_sim_data`: Selected columns of sim_data_filename including flux index and timepoints for training
	* `cols_sim_param`: Selected columns of sim_param_filename for training
	* `order_input_mids`: Order of CNN input MIDs as a list
	* `max_time`: Maximum timepoint to consider for training (hour)
	* `min_time`: Minimum timepoint to consider for training  (hour)
	* `delta_t_sim`: Time span between two timepoints in simulation (hour)
	* `n_metabs`: Number of metabolites used for training
	* `n_mids`: Number of isotopologues used for training
	* `metab_sources`: A list of reactions that produce a target metabolite
	* `target_flux`: A flux predicted by CNN, needs to be an element of metab_sources 
	* `name_rel_flux`: A name for the ratio of target_flux to metab_sources 
	* `non_mid_cols`: A list of input non-MID parameters, e.g., flux index and timepoints 
	* `hypothetical_MIDs`: If number of carbons in metabolites are different, for metabolites with fewer number of carbons, hypothetical_MIDs are defined as a list to make the input data consistent with CNN architecture.
	
	The information in the JSON file can also be found in the help function:
	```
	python run_metabolic_CNN.py --help
	```
	
1. Run metabolic CNN: The only input argument is a JSON file which is set to plasma serine uptake flux in glioma model configuration file `serine_plasma_glioma_config`, by default. To Run it with other config files:
	```
	python run_metabolic_CNN.py --config_file serine_denovo_glioma_config
	```
	Use the following config files to reproduce the results shown in the manuscript:
	* `serine_denovo_cortex_config`: train CNN to predict relative flux of glucose-derived serine synthesis in cortex
	* `serine_denovo_glioma_config`: train CNN to predict relative flux of glucose-derived serine synthesis in glioma
	* `serine_plasma_glioma_config`: train CNN to predict relative flux of plasma serine uptake in glioma
	* `gmp_denovo_glioma_config`: train CNN to predict relative flux of GMP _de novo_ synthesis in glioma
	
	This script performs the following steps:
	* Prepare simulated MIDs and fluxes for CNN model: This step creates an object of class `read_sim_data.py` and performs shuffling, splitting, normalizing, and reshaping simulated data described in help for all class functions.
	* Train CNN model on simulated data
	* Evaluate CNN performance on test simulated data: This step generates Fig. 4C, and 5C using `plotFunctions.py`.  
		
	You need to run `data_simulation` first to generate fluxes and MIDs and then place the generated files in `data/sim_data` folder. Alternatively, you can download simulated datasets from and move them to `data/sim_data` folder.
	
1. Predict _in vivo_ fluxes: Run `CNN_pred_mouse_patient.py` to predict relative fluxes in mice and patients using the trained model. 
	```
	python CNN_pred_mouse_patient.py --config_file pred_gmp_denovo_glioma_config
	```
	This script needs a JSON config file as an input argument which contains the following information:
    * `model_dir`: Directory of trained CNN model
    * `mice_dir`: Directory of mice MID data 
    * `mice_mid_name`: Column names of mice MID data (MIDs)
	* `mice_time`: Time of continuous infusion of uniformly labeled <sup>13</sup>C-glucose in mice
	* `patient_dir`: Directory of patient MID data 
	* `patient_mid_name`: Column names of mice MID data (MIDs)
	* `patient_time`: Directory of time of continuous infusion of uniformly labeled <sup>13</sup>C-glucose in mice
	* `order_input_mids`: Order of CNN input MIDs as a list
	* `hypothetical_MIDs`: If number of carbons in metabolites are different, for metabolites with fewer number of carbons, hypothetical_MIDs are defined as a list to make the input data consistent with CNN architecture.
	* `max_time`: Maximum timepoint to consider for training (hour)
	* `min_time`: Minimum timepoint to consider for training (hour)
	* `delta_t_sim`: Time span between two timepoints in simulation (hour)
	* `n_metabs`: Number of metabolites used for training
	* `n_mids`: Number of isotopologues used for training

	The information in the JSON file can be also found in the help function:
	```
	python CNN_pred_mouse_patient.py --help
	```
	Use the following config files to reproduce the results shown in the manuscript:
	* `pred_serine_denovo_cortex_config`: predict relative flux of glucose-derived serine synthesis in cortex for GBM12-, GBM38-, HF2303-bearing mice and patients
	* `pred_serine_denovo_glioma_config`: predict relative flux of glucose-derived serine synthesis in glioma for GBM12-, GBM38-, HF2303-bearing mice and patients
	* `pred_serine_plasma_glioma_config`: predict relative flux of plasma serine uptake in glioma for GBM12-, GBM38-, HF2303-bearing mice and patients
	* `pred_gmp_denovo_glioma_config`: predict relative flux of GMP _de novo_ synthesis in glioma for MMF-treated GBM38-bearing mice and patients
	* `pred_gmp_denovo_TRP_GBM38_config`: predict relative flux of GMP _de novo_ synthesis in TRP and GBM38 tumors