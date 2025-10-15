import json
from pathlib import Path
import os
import argparse

import tensorflow as tf
import tensorflow.keras as tfk
import pickle
import pandas as pd
import numpy as np

def main(args: argparse.Namespace):
    #%% Step 0: Set configurations
    
    # change working directory to source file location
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    
    # config filename without extension
    config_filename = args.config_file
    
    # directory of config file
    config_dir = "".join(['./', config_filename, '.json'])
    
    if not isinstance(config_dir, Path):
        config_dir = Path(config_dir)
    # Read config.json file
    with open(config_dir, 'r') as j:
        config_dict = json.loads(j.read())
        
    # check whether directory exists
    if not os.path.exists(config_dict['model_dir']):
        raise RuntimeError("CNN model output folder not found.")
    
    model_file = Path(config_dict['model_dir']) / 'cnn_model.h5'
    scaler_file = Path(config_dict['model_dir']) / 'scaler_attr.pkl'
    
    keep_col_data = config_dict['order_input_mids']
    
    n_t = int((config_dict['max_time'] - config_dict['min_time']) / 
            config_dict['delta_t_sim'] + 1)
    
    hypo_mids = config_dict['hypothetical_MIDs']
    
    with open(scaler_file, 'rb') as f:
        scaler = pickle.load(f)
    
    scaler_mean = scaler['scaler_mean']
    scaler_std = scaler['scaler_std']
    
    #%% Step 1: Load trained CNN model
    
    print('Physical GPUs:', tf.config.list_physical_devices('GPU'))
    
    model = tfk.models.load_model(model_file)
    
    #%% Step 2: Predict relative target flux in mice
    if config_dict['mice_dir'] != '':
        mice_files = os.listdir(config_dict['mice_dir'])
        mice_mid_names = pd.read_csv(config_dict['mice_mid_name'], 
                                    sep = '\t', header = None)
        
        mice_pred_dir = "".join([config_dict['model_dir'], '/',
                                    'mice', '_', config_filename])
        
        if not os.path.exists(mice_pred_dir):
            os.mkdir(mice_pred_dir)
            print("Folder %s created!" % mice_pred_dir)
        else:
            print("Folder %s already exists" % mice_pred_dir)
        
        for ps in range(len(mice_files)):
            mid_mc_dir = "".join([config_dict['mice_dir'], '/', mice_files[ps]])
            if mid_mc_dir[-3:] == 'csv':
                mid_mc = pd.read_csv(mid_mc_dir, sep = ',', header = 0, index_col = 0)
                mid_mc = mid_mc / 100
            elif mid_mc_dir[-3:] == 'txt':
                mid_mc = pd.read_csv(mid_mc_dir, sep = '\t', header = None)
                
            mid_mc.columns = mice_mid_names[0]
            mid_mc = mid_mc[keep_col_data]
            
            
            if hypo_mids != '':
                for i in range(len(hypo_mids)):
                    mid_mc[hypo_mids[i]] = 0
            
            # standardize mice data
            scaler_row_time = ((n_t - 1) / (config_dict['max_time']
                                            - config_dict['min_time']) * 
                            (config_dict['mice_time'] - config_dict['min_time']))
            scaler_row_time = int(scaler_row_time)
            
            scaled_mice_mids = ((mid_mc - scaler_mean[scaler_row_time, ]) / 
                                scaler_std[scaler_row_time, ])
            
            ML_in_mice_mids = np.reshape(scaled_mice_mids, (scaled_mice_mids.shape[0], 
                                                            config_dict['n_metabs'], 
                                                            config_dict['n_mids'], 
                                                            1))
                
            # model prediction
            ML_pred_mice_flux = model.predict(ML_in_mice_mids)
            
            pred_mice_df = pd.DataFrame(data = ML_pred_mice_flux, 
                                        columns = [config_filename])
                
            pred_mice_df.to_csv("".join([mice_pred_dir, '/', mice_files[ps][0:-4],
                                        '.txt']), sep = '\t')
    
    #%% Step 3: Predict relative target flux in patients
    if config_dict['patient_dir'] != '':
        patient_files = os.listdir(config_dict['patient_dir'])
        patient_mid_names = pd.read_csv(config_dict['patient_mid_name'], 
                                    sep = '\t', header = None)
        
        patient_time = pd.read_csv(config_dict['patient_time'], 
                                    sep = '\t', header = None)
        patient_time[1] = (patient_time[1] / 60).round(1)
        
        patient_pred_dir = "".join([config_dict['model_dir'], '/',
                                    'patients', '_', config_filename])
        
        if not os.path.exists(patient_pred_dir):
            os.mkdir(patient_pred_dir)
            print("Folder %s created!" % patient_pred_dir)
        else:
            print("Folder %s already exists" % patient_pred_dir)
        
        for ps in range(len(patient_files)):
            mid_mc_dir = "".join([config_dict['patient_dir'], '/', patient_files[ps]])
            if mid_mc_dir[-3:] == 'csv':
                mid_mc = pd.read_csv(mid_mc_dir, sep = ',', header = 0, index_col = 0)
                mid_mc = mid_mc / 100
            elif mid_mc_dir[-3:] == 'txt':
                mid_mc = pd.read_csv(mid_mc_dir, sep = '\t', header = None)
                
            mid_mc.columns = patient_mid_names[0]
            mid_mc = mid_mc[keep_col_data]
            
            if hypo_mids != '':
                for i in range(len(hypo_mids)):
                    mid_mc[hypo_mids[i]] = 0
            
            # standardize patient data
            scaler_row_time = ((n_t - 1) / (config_dict['max_time']
                                            - config_dict['min_time']) * 
                            (patient_time.loc[ps, 1] - config_dict['min_time']))
            scaler_row_time = int(scaler_row_time)
            
            scaled_patient_mids = ((mid_mc - scaler_mean[scaler_row_time, ]) / 
                                scaler_std[scaler_row_time, ])
            
            ML_in_patient_mids = np.reshape(scaled_patient_mids, 
                                            (scaled_patient_mids.shape[0], 
                                            config_dict['n_metabs'], 
                                            config_dict['n_mids'], 1))
                
            # model prediction
            ML_pred_patient_flux = model.predict(ML_in_patient_mids)
            
            pred_patient_df = pd.DataFrame(data = ML_pred_patient_flux, 
                                        columns = [config_filename])
                
            pred_patient_df.to_csv("".join([patient_pred_dir, '/', 
                                            patient_files[ps][0:3],
                                            '.txt']), sep = '\t')
    return

def parse_arguments(parser: argparse.ArgumentParser):

    parser.add_argument(
        "--config_file",
        type = str,
        default = "pred_gmp_denovo_TRP_GBM38_glioma_config",
        help = "Input arguments of the model are provided in a JSON file which contains the following information: \n"
        "model_dir: Directory of trained model \n"
	    "mice_dir:  Directory of mice MID data \n"
	    "mice_mid_name: Column names of mice MID data (MIDs) \n"
	    "mice_time: Time of continuous infusion of uniformly labeled 13C-glucose in mice \n"
        "patient_dir: Directory of patient MID data \n"
        "patient_mid_name: Column names of patient MID data \n"
        "patient_time: Time of continuous infusion of uniformly labeled 13C-glucose in patients \n"
	    "order_input_mids: Order of CNN input MIDs as a list \n"
	    "hypothetical_MIDs: If number of carbons in metabolites are different, for metabolites with fewer number of carbons, hypothetical_MIDs are defined as a list to make the input data consistent with CNN architecture. \n"
        "max_time: Maximum timepoint to consider for training (hour) \n"
	    "min_time: Minimum timepoint to consider for training  (hour) \n"
	    "delta_t_sim: time span between two timepoints in simulation (hour) \n"
	    "n_metabs: Number of metabolites used for training \n"
	    "n_mids: Number of isotopologues used for training \n"
        "Please see the default JSON file as an example."
    )
    
    arg = parser.parse_args()

    return arg

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Metabolic CNN - Prediction of relative _in vivo_ fluxes of metabolite sources",
        formatter_class=argparse.RawTextHelpFormatter
    )
    args = parse_arguments(parser)
    main(args)    
    