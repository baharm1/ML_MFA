import tensorflow as tf
import tensorflow.keras as tfk
from tensorflow.keras.optimizers import Adam
import datetime
import pickle
import pandas as pd
import numpy as np
from scipy.stats import linregress
import os
import argparse

from plotFunctions import plt_learning_curve, plt_target_pred, plot_regression
from read_sim_data import read_sim_data

def main(args: argparse.Namespace):

    #%% Step 0: Create output directory
    
    # config filename without extension
    config_filename = args.config_file
    
    # directory of config file
    config_dir = "".join(['./', config_filename, '.json'])
    
    # change working directory to source file location
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    
    # save time to include in output folder
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M")
    
    output_dir = "".join(['./', config_filename, '_', timestamp])
    
    # check whether directory already exists
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        print("Folder %s created!" % output_dir)
    else:
        print("Folder %s already exists" % output_dir)
    
    #%% Step 1: Prepare simulated fluxes and MIDs for CNN
    print('Preparing data for CNN')
    # create simulated data object
    data_obj = read_sim_data(config_dir)
    
    # split data to train and test dataframes
    train_data, test_data = data_obj.split_sim_data(data_obj.combined_data)
    
    # fit normalization to training features
    scaler = data_obj.fit_scaler(train_data)
    
    # normalize features 
    train_data = data_obj.transform_scaler_sim_data(scaler, train_data)
    test_data = data_obj.transform_scaler_sim_data(scaler, test_data)
    
    # shuffle and split features and target
    train_x, train_y = data_obj.split_features_target(train_data)
    x_test, y_test = data_obj.split_features_target(test_data)
    
    # split train and validation
    x_train, x_val = data_obj.split_sim_data(train_data, split_frac = 0.2)
    x_train, y_train = data_obj.split_features_target(x_train)
    x_val, y_val = data_obj.split_features_target(x_val)
    
    # reshape mids for cnn
    x_train = data_obj.reshape_features(x_train)
    x_val = data_obj.reshape_features(x_val) 
    x_test = data_obj.reshape_features(x_test)
    
    val_data = (x_val, y_val)
    
    # save scaler attributes
    scaler_mean = np.reshape(scaler.mean_, 
                            (data_obj.n_t, data_obj.n_metabs * data_obj.n_mids))
    scaler_std = np.reshape(np.sqrt(scaler.var_), 
                            (data_obj.n_t, data_obj.n_metabs * data_obj.n_mids))
    
    scaler_mean[scaler_std == 0] = 0
    scaler_std[scaler_std == 0] = 1
    
    scaler_attr = {"scaler_mean": scaler_mean, "scaler_std": scaler_std}
    
    with open("".join([output_dir, '/scaler_attr.pkl']), 'wb') as f:
        pickle.dump(scaler_attr, f)
        
    #%% Step 2: Train a CNN Model with one conv2d and one conv1d followed by a fully connected neural net
    
    print('Physical GPUs:', tf.config.list_physical_devices('GPU'))
    print('Start training')
    kernel_initializer = tfk.initializers.HeUniform(seed = 123)
    
    input_layer = tfk.Input(shape = (data_obj.n_metabs, data_obj.n_mids, 1))
    
    conv2d1 = tfk.layers.Conv2D(filters = 24, 
                                kernel_size = (data_obj.n_metabs, 1), 
                                activation = 'relu', 
                                kernel_initializer = kernel_initializer,
                                use_bias = False)(input_layer)
    
    batch1 = tfk.layers.BatchNormalization()(conv2d1)
    
    conv1d1 = tfk.layers.Conv1D(filters = 40, 
                                kernel_size = (data_obj.n_mids, ), 
                                activation = 'relu', 
                                kernel_initializer = kernel_initializer,
                                use_bias = False)(batch1)
    
    batch2 = tfk.layers.BatchNormalization()(conv1d1)
    
    flat1 = tfk.layers.Flatten()(batch2)
    
    dense1 = tfk.layers.Dense(32, 
                            activation = 'relu', 
                            kernel_initializer = kernel_initializer,
                            use_bias = False)(flat1)
    
    batch3 = tfk.layers.BatchNormalization()(dense1)
    
    dense2 = tfk.layers.Dense(16, 
                            activation = 'relu', 
                            kernel_initializer = kernel_initializer,
                            use_bias = False)(batch3)
    
    batch4 = tfk.layers.BatchNormalization()(dense2)
    
    output_layer = tfk.layers.Dense(1, 
                                    kernel_initializer = kernel_initializer, 
                                    activation = 'sigmoid')(batch4)
    
    model = tfk.Model(inputs = input_layer, outputs = output_layer)
    
    model.compile(loss = 'mse',  
                optimizer = Adam(learning_rate = 0.004, epsilon = 0.4), 
                metrics = ['mse'])
    
    print(model.summary())
    
    tfk.utils.plot_model(model, "".join([output_dir, '/cnn_architecture.png']), 
                        show_shapes = True, show_layer_names = True)
    
    es = []
    es.append(tfk.callbacks.EarlyStopping(monitor = 'val_loss', 
                                        patience = 8, 
                                        min_delta = 1e-6))
    es.append(tfk.callbacks.ReduceLROnPlateau(monitor = 'val_loss', 
                                            factor = 0.9,
                                            patience = 5, 
                                            min_lr = 0.001))                                 
    
    history = model.fit(x_train, y_train, 
                        epochs = 1000, 
                        batch_size = 256,
                        validation_data = val_data, 
                        callbacks = es)
    
    plt_learning_curve(history, 
                    filename = "".join([output_dir, '/learning_curve']))
    
    model.save("".join([output_dir, '/cnn_model.h5']))
    
    print("CNN model training finished.")
    with open("".join([output_dir, '/history.pkl']), 'wb') as f:
        pickle.dump(history, f)
    
    #%% Step 3: Evaluate CNN model on train and test data
    print('Visualizing CNN evaluation')
    pred_train = np.reshape(model.predict(x_train), newshape = y_train.shape)
    pred_test = np.reshape(model.predict(x_test), newshape = y_test.shape)
    
    plt_target_pred(y_train, pred_train, istrain = True,
                    filename = "".join([output_dir, '/corr_plot_train']))
    plt_target_pred(y_test, pred_test, istrain = False,
                    filename = "".join([output_dir, '/corr_plot_test']))
    
    pred_corr = linregress(y_train, pred_train)
    pred_corr = {'pearsonr': pred_corr.rvalue, 
                'slope': pred_corr.slope, 'intercept': pred_corr.intercept, 
                'slope_stderr': pred_corr.stderr, 
                'intercept_stderr': pred_corr.intercept_stderr,
                'R2': pred_corr.rvalue ** 2, 'pval': pred_corr.pvalue}
    pred_corr = pd.DataFrame([pred_corr])
    pred_corr.to_csv("".join([output_dir, '/corr_train.csv']))
    
    pred_corr = linregress(y_test, pred_test)
    pred_corr = {'pearsonr': pred_corr.rvalue, 
                'slope': pred_corr.slope, 'intercept': pred_corr.intercept, 
                'slope_stderr': pred_corr.stderr, 
                'intercept_stderr': pred_corr.intercept_stderr,
                'R2': pred_corr.rvalue ** 2, 'pval': pred_corr.pvalue}
    pred_corr = pd.DataFrame([pred_corr])
    pred_corr.to_csv("".join([output_dir, '/corr_test.csv']))
    
    plot_regression(y_test, pred_test,
                    filename = "".join([output_dir, 
                                        '/density_plot_test']))
    return

def parse_arguments(parser: argparse.ArgumentParser):

    parser.add_argument(
        "--config_file",
        type = str,
        default = "serine_plasma_glioma_config",
        help = "Input arguments of the model are provided in a JSON file which contains the following information: \n"
        "sim_data_dir: Directory of simulated data that contains sim_data_filename and sim_param_filename \n"
	    "sim_data_filename: Filename of simulated data that contains MIDs of balanced metabolites \n"
	    "sim_param_filename: Filename of simulated parameters that contains MIDs of input metabolites and fluxes \n"
	    "cols_sim_data: Selected columns of sim_data_filename including flux index and timepoints for training \n"
	    "cols_sim_param: Selected columns of sim_param_filename for training \n"
	    "order_input_mids: Order of MIDs as inputs of CNN \n"
	    "max_time: Maximum timepoint to consider for training (hour) \n"
	    "min_time: Minimum timepoint to consider for training  (hour) \n"
	    "delta_t_sim: time span between two timepoints in simulation (hour) \n"
	    "n_metabs: Number of metabolites used for training \n"
	    "n_mids: Number of isotopologues used for training \n"
        "metab_sources: A list of reactions that produce a target metabolite \n"
        "target_flux: A flux predicted by CNN, needs to be an element of metab_sources \n"
        "name_rel_flux: A name for the ratio of target_flux to metab_sources \n"
        "non_mid_cols: A list of input non-MID parameters, e.g., flux index and timepoints \n"
        "hypothetical_MIDs: If number of carbons in metabolites are different, for metabolites with fewer number of carbons, hypothetical_MIDs are defined as a list to make the input data consistent with CNN architecture. \n"
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
    