import tensorflow as tf
import tensorflow.keras as tfk
from tensorflow.keras.optimizers import Adam
import datetime
import pickle
import pandas as pd
import numpy as np
from scipy.stats import linregress
import os

from plotFunctions import plt_learning_curve, plt_target_pred, plot_regression
from read_sim_data import read_sim_data

#%% Step 0: Create output directory

# config filename without extension
config_filename = 'serine_plasma_glioma_config'

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


with open("".join([output_dir, '/history.pkl']), 'wb') as f:
    pickle.dump(history, f)

#%% Step 3: Evaluate CNN model on train and test data

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
