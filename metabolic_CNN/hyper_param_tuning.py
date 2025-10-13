import tensorflow.keras as tfk
from tensorflow.keras.optimizers import Adam
import matplotlib.pyplot as plt
import time
import os
import datetime
import optuna
from joblib import parallel_backend
from sklearn.metrics import r2_score
import plotly.io as io

from read_sim_data import read_sim_data
#%% set seed
import random
random.seed(123)

#%% Step 0: Create output directory

# config filename without extension
config_filename = 'gmp_denovo_glioma_config'

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
#%%
no_mids = 5
no_metabs = 4

#%%

def cv_optuna(trial):
    neurons = trial.suggest_int('neurons', low = 8, high = 12, step = 4)
    
    layers = trial.suggest_int('layers', low = 1, high = 3, step = 1)
    
    learning_rate = trial.suggest_float('learning_rate', low = 0.001, 
                                        high = 0.01, log = True)
    epsilon = trial.suggest_float('epsilon', low = 0.001, high = 1, log = True)
    
    bias = trial.suggest_float('bias', low = 0, high = 0.1, step = 0.05)    
    bias_initializer = tfk.initializers.Constant(value = bias)
    
    reg = trial.suggest_float('reg', low = 0.001, high = 0.1, log = True)
    kernel_regularizer = tfk.regularizers.L2(l2 = reg)
    
    filters = trial.suggest_int('filters', low = 4, high = 16, step = 2)
    
    dropout_rate = trial.suggest_float('dropout_rate', low = 0.1, 
                                       high = 0.3, step = 0.1)
    
    batch_size = trial.suggest_int('batch_size', low = 64, high = 512, step = 64)

    model = tfk.Sequential()
    model.add(tfk.layers.Conv2D(filters, (data_obj.n_metabs, 1), 
                                activation = 'relu', 
                                kernel_initializer = 'he_uniform', 
                                kernel_regularizer = kernel_regularizer,
                                bias_initializer = bias_initializer,
                                input_shape = (data_obj.n_metabs, 
                                               data_obj.n_mids, 1)))
    model.add(tfk.layers.Flatten())
    model.add(tfk.layers.BatchNormalization())
    model.add(tfk.layers.Dropout(dropout_rate))
    for _ in range(layers):
        model.add(tfk.layers.Dense(neurons, activation = 'relu', 
                                   kernel_initializer = 'he_uniform',
                                   kernel_regularizer = tfk.regularizers.L2(l2 = reg),
                                   bias_initializer = bias_initializer))
        model.add(tfk.layers.BatchNormalization())        
        model.add(tfk.layers.Dropout(dropout_rate))
        
    model.add(tfk.layers.Dense(1, kernel_initializer = 'he_uniform', 
                               kernel_regularizer = tfk.regularizers.L2(l2 = reg)),
                               activation = 'sigmoid')
                                   
    model.compile(loss = 'mean_squared_error', 
                  optimizer = Adam(learning_rate = learning_rate, 
                                   epsilon = epsilon), 
                  metrics = ['mse'])
    print(model.summary())
        
    
    es = tfk.callbacks.EarlyStopping(monitor='val_loss', patience = 5, 
                                     start_from_epoch = 50)
    
    with parallel_backend(backend='threading', n_jobs = -1):
        history = model.fit(x_train, y_train, epochs = 200, 
                            batch_size = batch_size,
                            validation_data = val_data, callbacks = [es])
    
    # plot learning curves
    plt.xlabel('Epoch')
    plt.ylabel('MSE')
    plt.plot(history.history['loss'], label = 'train')
    plt.plot(history.history['val_loss'], label = 'val')
    plt.legend()
    plt.show()
    
    pred = model.predict(x_val)
    coefficient_of_dermination = r2_score(y_val, pred)

    return coefficient_of_dermination

nn_bo = optuna.create_study(sampler = optuna.samplers.TPESampler(), 
                            direction = 'maximize')
start = time.time()
nn_bo.optimize(cv_optuna, n_trials = 10)
print('It takes %s minutes' % ((time.time() - start)/60))

#%%
import pickle
with open('optuna_cnn304min.pkl', 'wb') as f:
    pickle.dump(nn_bo, f)
#%%
fig = optuna.visualization.plot_optimization_history(nn_bo)
fig.write_image("optuna_cnn.pdf")
#%%
fig = optuna.visualization.plot_parallel_coordinate(nn_bo)
fig.write_image('optuna_cnn_parallel_coordinate.pdf')
io.renderers.default = 'browser'
fig.show()
#%%
fig = optuna.visualization.plot_contour(nn_bo)
fig.write_image('optuna_cnn_plot_contour.pdf')
#%%
fig = optuna.visualization.plot_contour(nn_bo, params = ["layers", "neurons"])
fig.write_image('optuna_cnn_contour_neurons_layers.pdf')
#%%
fig = optuna.visualization.plot_contour(nn_bo, params = ["learning_rate", "epsilon"])
fig.write_image('optuna_cnn_contour_lr_epsilon.pdf')
#%%
fig = optuna.visualization.plot_contour(nn_bo, params = ["batch_size", "epochs"])
fig.write_image('optuna_cnn_contour_batch_size_epochs.pdf')
#%%
fig = optuna.visualization.plot_contour(nn_bo, params = ["reg", "epsilon"])
fig.write_image('optuna_cnn_contour_reg_epsilon.pdf')
#%%
fig = optuna.visualization.plot_contour(nn_bo, params = ["activation", "reg"])
fig.write_image('optuna_cnn_contour_activation_reg.pdf')
#%%
fig = optuna.visualization.plot_slice(nn_bo)
fig.write_image('optuna_cnn_plot_slice.pdf')
#%%
fig = optuna.visualization.plot_param_importances(nn_bo)
fig.write_image('optuna_cnn_plot_param_importance.pdf')
io.renderers.default = 'browser'
fig.show()
