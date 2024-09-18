# https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
import pandas as pd
import numpy as np
import tensorflow.keras as tfk
from tensorflow.keras.optimizers import Adam
import matplotlib.pyplot as plt
import time
import optuna
from joblib import parallel_backend
from sklearn.model_selection import GroupKFold, cross_val_score
from scikeras.wrappers import KerasRegressor
from sklearn.metrics import r2_score
from sklearn.preprocessing import MinMaxScaler
import plotly.io as io
#%% set seed
import random
random.seed(123)
#%%
# Read the data file
data = pd.read_csv('..//Data Simulation//simulated_data_v5.csv')
# print(data.columns) # view the columns
#%%
# Keep the relevant features
columns_to_keep = ['index', 'time',
                   'IMP1', 'IMP2', 'IMP3', 'IMP4', 'IMP5',
                   'GMP1', 'GMP2', 'GMP3', 'GMP4', 'GMP5',
                   'R5P1', 'R5P2', 'R5P3', 'R5P4', 'R5P5',
                   'GDP1', 'GDP2', 'GDP3', 'GDP4', 'GDP5', 
                   'ratio']
data = data[columns_to_keep]
data = data.loc[data['time'] >= 2]
data = data.drop(['time'], axis = 1)

#%%
data = data.loc[data['time'] == 4]
data = data.drop(['time'], axis = 1)

#%%
# Divide into test and train sets such that the test set had unique conditions
def split_train_test_data(data, test_size = 0.15):
    # Divide into test and train sets such that the test set had unique conditions
    index = data['index'].unique() # simulations with the same fluxes      
    test_index = random.sample(list(index), round(len(index) * test_size))
    #def sample_function():
    #    return 0.5
    #index = random.shuffle(list(index), sample_function)
    #test_index = index[0:round(len(index) * test_size)]
    test_data = data.loc[data['index'].isin(test_index)]
    train_data = data.loc[~(data['index'].isin(test_index))] # including train and validation
    
    return(train_data, test_data)

def split_data_features_outcome(data):
    # shuffle the DataFrame rows
    data = data.sample(frac = 1)

    # Divide the data into features and outcome
    features = data.drop(columns = ['index', 'ratio'])
    outcome = data['ratio'].to_numpy()
    flux_index = data['index'].to_numpy()
    features = features.to_numpy()
    
    return(features, outcome, flux_index)

def scale_min_max(data):
    scaler = MinMaxScaler()
    scaler.fit(data)
    data = scaler.transform(data)
    return data

train_data, test_data = split_train_test_data(data)
train_x, train_y, train_ind = split_data_features_outcome(train_data)
x_test, y_test, ind_test = split_data_features_outcome(test_data)

x_train, x_val = split_train_test_data(train_data, test_size = 0.15)
x_train, y_train, ind_train = split_data_features_outcome(x_train)
x_val, y_val, ind_val = split_data_features_outcome(x_val)
#%%
train_x = np.transpose(np.transpose(train_x) / np.max(train_x, axis = 1))
x_train = np.transpose(np.transpose(x_train) / np.max(x_train, axis = 1))
x_test = np.transpose(np.transpose(x_test) / np.max(x_test, axis = 1))
x_val = np.transpose(np.transpose(x_val) / np.max(x_val, axis = 1))

#%%
no_mids = 5
no_metabs = 4

def my_scale_max(data):
    data = np.reshape(data, (data.shape[0], no_metabs, no_mids, 1))
    data = data / np.max(data, axis = 1, keepdims = True)
    #for i in range(data.shape[0]):
       # maxim = np.max(data[i, :, :, 0], axis = 0)
      #  data[i, :, :, 0] = data[i, :, :, 0] / maxim
     #   minim = np.min(data[i, :, :, 0], axis = 0)
        
      #  data[i, :, :, 0] = (data[i, :, :, 0] - minim) / (maxim - minim)
    return data
    
train_x = my_scale_max(train_x)
x_train = my_scale_max(x_train)
x_val = my_scale_max(x_val)
x_test = my_scale_max(x_test)

val_data = (x_val, y_val)
#%%
train_x = scale_min_max(train_x)
x_train = scale_min_max(x_train)
x_val = scale_min_max(x_val)
x_test = scale_min_max(x_test)

#%%


for i in range(x_train.shape[0]):
    v = np.max(x_train[i,:], axis = 0)
    x_train[i,:] = x_train[i,:]/v
    
for i in range(x_val.shape[0]):
    v = np.max(x_val[i,:], axis = 0)
    x_val[i,:] = x_val[i,:]/v

for i in range(x_test.shape[0]):
    v = np.max(x_test[i,:], axis = 0)
    x_test[i,:] = x_test[i,:]/v    

for i in range(train_x.shape[0]):
    v = np.max(train_x[i,:], axis = 0)
    train_x[i,:] = train_x[i,:]/v
val_data = (x_val, y_val)


#%%

train_x = np.reshape(train_x, (train_x.shape[0], no_metabs * no_mids))
x_train = np.reshape(x_train, (x_train.shape[0], no_metabs * no_mids))
x_val = np.reshape(x_val, (x_val.shape[0], no_metabs * no_mids))
x_test = np.reshape(x_test, (x_test.shape[0], no_metabs * no_mids))
val_data = (x_val, y_val)
#%%

def cv_optuna(trial):
    neurons = trial.suggest_int('neurons', low = 8, high = 12, step = 4)
    
    layers = trial.suggest_int('layers', low = 1, high = 3, step = 1)
    
    learning_rate = trial.suggest_float('learning_rate', low = 0.001, high = 0.01, log = True)
    epsilon = trial.suggest_float('epsilon', low = 0.001, high = 1, log = True)
    
    bias = trial.suggest_float('bias', low = 0, high = 0.1, step = 0.05)    
    bias_initializer = tfk.initializers.Constant(value = bias)
    
    reg = trial.suggest_float('reg', low = 0.001, high = 0.1, log = True)
    kernel_regularizer = tfk.regularizers.L2(l2 = reg)
    
    filters = trial.suggest_int('filters', low = 4, high = 16, step = 2)
    
    dropout_rate = trial.suggest_float('dropout_rate', low = 0.1, high = 0.3, step = 0.1)
    
    batch_size = trial.suggest_int('batch_size', low = 64, high = 512, step = 64)

    model = tfk.Sequential()
    model.add(tfk.layers.Conv2D(filters, (no_metabs, 1), 
                                activation = 'relu', 
                                kernel_initializer = 'he_uniform', 
                                kernel_regularizer = kernel_regularizer,
                                bias_initializer = bias_initializer,
                                input_shape = (no_metabs, no_mids, 1)))
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
        
    model.add(tfk.layers.Dense(1, kernel_initializer = 'glorot_uniform', 
                               kernel_regularizer = tfk.regularizers.L2(l2 = reg)))
                               #activation = 'sigmoid')) 
                                   
    model.compile(loss = 'mean_squared_error', 
                  optimizer = Adam(learning_rate = learning_rate, epsilon = epsilon), 
                  metrics = ['mse'])
    print(model.summary())
        
    
    es = tfk.callbacks.EarlyStopping(monitor='val_loss', patience = 5, 
                                     start_from_epoch = 50), 
                                     #restore_best_weights = True)
    #val_data = (x_val, y_val)
    with parallel_backend(backend='threading', n_jobs = -1):
        history = model.fit(x_train, y_train, epochs = 200, 
                            batch_size = batch_size,
                            validation_data = val_data, callbacks = [es])
    
    # plot learning curves
    #plt.title('Learning Curves')
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
io.renderers.default = 'browser'
fig.show()
#fig.write_image('optuna_cnn_parallel_coordinate.pdf')
#%%
fig = optuna.visualization.plot_contour(nn_bo)#, params = ["layers", "neurons"])
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

#%%
y_train = y_train * 100
y_val = y_val * 100
val_data = (x_val, y_val)
y_test = y_test * 100
#%%
bias_initializer = tfk.initializers.Constant(value = 0.001)
kernel_regularizer = tfk.regularizers.L2(l2 = 0.001)
kernel_initializer = tfk.initializers.HeUniform(seed = 123)
model = tfk.Sequential()
model.add(tfk.layers.Conv2D(12, (no_metabs, 1), activation = 'relu', 
                            kernel_initializer = kernel_initializer, 
                            kernel_regularizer = kernel_regularizer, use_bias = False,
                            #bias_initializer = bias_initializer,
                            input_shape = (no_metabs, no_mids, 1)))

model.add(tfk.layers.BatchNormalization())
model.add(tfk.layers.Conv1D(12, (no_mids,), activation = 'relu', 
                            kernel_initializer = kernel_initializer, 
                            kernel_regularizer = kernel_regularizer, use_bias = False))
                            
model.add(tfk.layers.Flatten())
model.add(tfk.layers.BatchNormalization())
model.add(tfk.layers.Dropout(0.1))

model.add(tfk.layers.Dense(30, activation = 'relu', 
                           kernel_initializer = kernel_initializer,
                           kernel_regularizer = kernel_regularizer, use_bias = False))#,
#bias_initializer = bias_initializer))
model.add(tfk.layers.BatchNormalization())        
model.add(tfk.layers.Dropout(0.1))

model.add(tfk.layers.Dense(16, activation = 'relu', 
                           kernel_initializer = kernel_initializer,
                           kernel_regularizer = kernel_regularizer, use_bias = False))#,
#bias_initializer = bias_initializer))
model.add(tfk.layers.BatchNormalization())        
model.add(tfk.layers.Dropout(0.1))

model.add(tfk.layers.Dense(1, kernel_initializer = kernel_initializer, 
                           kernel_regularizer = kernel_regularizer,
                           activation = 'hard_sigmoid'))
model.compile(loss = 'mean_squared_error', 
              optimizer = Adam(learning_rate = 0.004, epsilon = 0.4), 
              metrics = ['mse'])
print(model.summary())

es = []
es.append(tfk.callbacks.EarlyStopping(monitor = 'val_loss', patience = 5, restore_best_weights = True))#, start_from_epoch = 50)) 
                                 #restore_best_weights = True)
#es.append(tfk.callbacks.ReduceLROnPlateau(monitor = 'val_loss', factor = 0.9,
#                              patience = 4, min_lr = 0.001))                                 
#val_data = (x_val, y_val)
with parallel_backend(backend = 'threading', n_jobs = -1):
    history = model.fit(x_train, y_train, epochs = 150, 
                        batch_size = 128,
                        validation_data = val_data, callbacks = es)
#%%
# save model to file
#model.save('model_hard_sigmoid.h5')
# plot learning curves
plt.title('Learning Curves')
plt.xlabel('Epoch')
plt.ylabel('Loss (MSE)')
plt.plot(history.history['loss'], label = 'Train')
plt.plot(history.history['val_loss'], label = 'Validation')
plt.legend()
#plt.savefig('learning_curve_cnn_hard_sigmoid.pdf', format = 'pdf', bbox_inches = 'tight')
plt.show()

#%%
# load the model from file
model = tfk.models.load_model('modelcnn.h5')

#%%


#%% model visualization
import visualkeras

from tensorflow.keras.layers import Dense, Conv2D, Flatten, Dropout, Conv1D, BatchNormalization
from collections import defaultdict

color_map = defaultdict(dict)
color_map[Conv2D]['fill'] = '#AE3C60'
color_map[BatchNormalization]['fill'] = '#F3C33C'
color_map[Dropout]['fill'] = 'darkcyan'
color_map[Conv1D]['fill'] = 'turquoise'
color_map[Dense]['fill'] = 'skyblue'
color_map[Flatten]['fill'] = '#255E79'

visualkeras.layered_view(model, legend=True, to_file = 'model_2.png', scale_xy=20, scale_z=1, color_map=color_map, spacing = 20)

#%%
from keras.utils.vis_utils import plot_model
plot_model(model, to_file='model_plot.png', show_shapes = True, show_layer_names = True)

#%%
fig, ax = plt.subplots(figsize = (5, 5))
pred = model.predict(x_test)
plt.scatter(y_test, pred)
b, a = np.polyfit(y_test, pred, deg = 1)

# Create sequence of 100 numbers from 0 to 100 
xseq = np.linspace(0, 1, num = 100)
# Plot regression line
ax.plot(xseq, a + b * xseq, color = "k", lw = 2.5);
coefficient_of_dermination = r2_score(y_test, pred)
plt.title("Test R2 = %0.2f " % coefficient_of_dermination)
plt.xlabel('actual')
plt.ylabel('pred')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.show()
#%%
pred = model.predict(x_train)
fig, ax = plt.subplots(figsize = (5, 5))
plt.scatter(y_train, pred)
b, a = np.polyfit(y_train, pred, deg = 1)

# Create sequence of 100 numbers from 0 to 100 
xseq = np.linspace(0, 1, num = 100)
# Plot regression line
ax.plot(xseq, a + b * xseq, color = "k", lw = 2.5);
coefficient_of_dermination = r2_score(y_train, pred)
plt.title("Train R2 = %0.2f " % coefficient_of_dermination)
plt.xlabel('actual')
plt.ylabel('pred')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.show()

#%%
data_mice = pd.read_csv('mice_data_MMF.csv')

# Keep the relevant features
columns_to_keep = ['IMP1', 'IMP2', 'IMP3', 'IMP4', 'IMP5',
                   'GMP1', 'GMP2', 'GMP3', 'GMP4', 'GMP5',
                   'R5P1', 'R5P2', 'R5P3', 'R5P4', 'R5P5',
                   'GDP1', 'GDP2', 'GDP3', 'GDP4', 'GDP5']
data_mice = data_mice[columns_to_keep].to_numpy()

data_mice = np.reshape(data_mice, (data_mice.shape[0], no_metabs, no_mids, 1))
for i in range(data_mice.shape[0]):
    maxim = np.max(data_mice[i, :, :, 0], axis = 0)
    maxim[maxim == 0] = 1e-15
    
    data_mice[i, :, :, 0] = data_mice[i, :, :, 0] / maxim

pred_mice = model.predict(data_mice)
np.savetxt("pred_mice_hard_sigmoid.csv", pred_mice, delimiter=",")
#%%
data_human = pd.read_csv('patient_unavg.csv')
#data_human = pd.read_csv('human_data.csv')
# Keep the relevant features
columns_to_keep = ['IMP1', 'IMP2', 'IMP3', 'IMP4', 'IMP5',
                   'GMP1', 'GMP2', 'GMP3', 'GMP4', 'GMP5',
                   'R5P1', 'R5P2', 'R5P3', 'R5P4', 'R5P5',
                   'GDP1', 'GDP2', 'GDP3', 'GDP4', 'GDP5']
data_human = data_human[columns_to_keep].to_numpy()

data_human = np.reshape(data_human, (data_human.shape[0], no_metabs, no_mids, 1))
for i in range(data_human.shape[0]):
    maxim = np.max(data_human[i, :, :, 0], axis = 0)
    maxim[maxim == 0] = 1e-15
    data_human[i, :, :, 0] = data_human[i, :, :, 0] / maxim

pred_human = model.predict(data_human)
np.savetxt("pred_human_unavg_hard_sigmoid.csv", pred_human, delimiter=",")
#%%
import scipy
def plot_regression(x, y, title):
    '''
    Plot the predicted vs ground truth figure. Include R2 and fit line.
    '''
    # Calculate the point density
    xy = np.vstack((x,y))
    z = scipy.stats.gaussian_kde(xy)(xy)

    fig, ax = plt.subplots()
    den = ax.scatter(x, y, c = z, s = 1)

    # regression
    xseq = np.linspace(x.min(), x.max(), num = 100)
    
    b, a = np.polyfit(x, y, deg = 1)
    ax.plot(xseq, a + b * xseq, color = "k", lw = 1.5)
    #ax.plot(xseq, xseq, color = "k", lw = 1, linestyle = 'dashed')
    #res = scipy.stats.linregress(x, y)
    #plt.text(0.08, 0.82, f'y = {res.slope:.3f}x + {res.intercept:.3f}', fontsize = 12, transform = ax.transAxes)
    #plt.text(0.6, 0.20, f'R-squared: {res.rvalue:.2f}' , fontsize = 11, transform = ax.transAxes)
    #plt.text(0.6, 0.10, f'n={len(x)}', fontsize = 11, transform = ax.transAxes)
    #plt.title(title, fontsize=13)
    coefficient_of_dermination = r2_score(x, y)
    plt.title(title + " (R2 = %0.2f)" % coefficient_of_dermination, fontsize = 13)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    ax.set_xlabel('Target', fontsize=13)
    ax.set_ylabel('Predicted', fontsize=13)

    fig.patch.set_facecolor('white')
    fig.colorbar(den, label = 'Density')
    #plt.savefig(title + '.pdf', format = 'pdf', bbox_inches = 'tight')

    #print (f'Printed plot {title=} (n={len(x)}) with R-squared: {res.rvalue:.2f}, slope: {res.slope:.3f}, and intercept: {res.intercept:.3f}')

pred_train = np.reshape(model.predict(x_train), newshape=y_train.shape)
plot_regression(y_train, pred_train, 'Training Dataset')
pred_val = np.reshape(model.predict(x_val), newshape=y_val.shape)
plot_regression(y_val, pred_val, 'Validation Dataset')
pred_test = np.reshape(model.predict(x_test), newshape=y_test.shape)
plot_regression(y_test, pred_test, 'Test Dataset')

#%% 
model.evaluate(x_val, y_val)
model.evaluate(x_test, y_test)
#%%
outputs = [model.layers[i].output for i in range(8)]

newModel = tfk.models.Model(model.inputs, outputs)
sample = x_train[1, :, :, :]

sample = np.reshape(sample, (1, 4, 5, 1))

newModel.predict(sample)

# save model to file
model.save('model.h5')


#%%
                                
history = model.fit(x_train, y_train, epochs = 200, batch_size = 128, 
                        validation_data = val_data, callbacks = [es])

# plot learning curves
plt.title('Learning Curves')
plt.xlabel('Epoch')
plt.ylabel('MSE')
plt.plot(history.history['loss'], label='train')
plt.plot(history.history['val_loss'], label='val')
plt.legend()
plt.show()
#%%

plt.scatter(y_test, model.predict(x_test))
plt.show()

plt.scatter()

#%%
def cv_optuna(trial):
    neurons = trial.suggest_int('neurons', low = min(neurons), 
                                high = max(neurons), 
                                step = min(neurons)/2)
    layers = trial.suggest_int('layers', low = min(layers), 
                               high = max(layers), step = 1)
    optimizer = trial.suggest_categorical('optimizer', optimizer)
    #learning_rate = trial.suggest_float('learning_rate', 
    #                                    low = min(learning_rate),
    #                                    high = max(learning_rate), 
    #                                    log = True)
    #optimizer = set_optimizer_lr(learning_rate, optimizer)
    activation = trial.suggest_categorical('activation', activation)
    kernel_initializer = trial.suggest_categorical('kernel_initializer', 
                                                   kernel_initializer)
    batch_size = trial.suggest_int('batch_size', low = min(batch_size), 
                                   high = max(batch_size), 
                                   step = min(batch_size))
    epochs = trial.suggest_int('epochs', low = min(epochs), 
                               high = max(epochs), 
                               step = min(epochs))
    
    def create_model():
        model = tfk.Sequential()
        # first hidden layer
        model.add(tfk.layers.Dense(neurons, activation = activation, 
                                   kernel_initializer = kernel_initializer, 
                                   input_shape = (nfeatures, )))
        # additional hidden layers
        for _ in range(layers - 2):
            model.add(tfk.layers.Dense(neurons, activation = activation,
                                       kernel_initializer = kernel_initializer))
        # output layer    
        model.add(tfk.layers.Dense(1, kernel_initializer = kernel_initializer))
        
        # summarize the model
        #print(model.summary())
        model.compile(loss = 'mean_squared_error', optimizer = optimizer, 
                      metrics = ['mse'])
        return model
    
    estimator = KerasRegressor(model = create_model, 
                               batch_size = batch_size, 
                               epochs = epochs, verbose = 0)
    cv_gfk = GroupKFold(n_splits = 5) # 5-fold split for cross validation   
    with parallel_backend(backend='threading', n_jobs=-1):
        cv_score = cross_val_score(estimator, train_x, train_y, 
                                   scoring = 'neg_mean_squared_error', 
                                   n_jobs = -1, 
                                   groups = train_ind, cv = cv_gfk).mean()
    return cv_score

def bayes_optuna(self):
    
    nn_bo = optuna.create_study(sampler = optuna.samplers.TPESampler(), 
                                direction = 'maximize')
    start = time.time()
    nn_bo.optimize(cv_optuna, n_trials = 100, timeout = 5400)
    print('It takes %s minutes' % ((time.time() - start)/60))
    return nn_bo
#%%
# Benchmark the task with a linear SVR model
# The M+0 is removed to remove the linear relationship between features
#from sklearn.svm import LinearSVR
#from sklearn.metrics import mean_squared_error
#linear_model = LinearSVR(C = 100, tol = 0.005, max_iter = 2000)
#linear_model.fit(np.reshape(X_train[:,1:,:],(-1,20),'F'),y_train)
#y_train_pred_linear = linear_model.predict(np.reshape(X_train[:,1:,:],(-1,20),'F'))
#print(mean_squared_error(y_train, y_train_pred_linear))


import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten, Conv2D, BatchNormalization, Dropout

kernel_init = tf.keras.initializers.HeNormal()
bias_init = tf.keras.initializers.Constant(value=0.1)
reg = tf.keras.regularizers.L2(l2=0.001)
opt = tf.keras.optimizers.Adam(learning_rate=0.005, epsilon = 0.1)
callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=2, min_delta = 0.002)
drop = 0.3
model = Sequential([
    Conv2D(8, (no_metabs, 1), activation='relu', input_shape=(no_metabs, no_mids, 1),
           kernel_initializer=kernel_init, bias_initializer=bias_init,
           kernel_regularizer=reg),
    BatchNormalization(),
    Flatten(),
    Dense(32, activation='relu', kernel_initializer=kernel_init, bias_initializer=bias_init,kernel_regularizer=reg),
    BatchNormalization(),
    Dropout(drop),
    Dense(32, activation='relu', kernel_initializer=kernel_init, bias_initializer=bias_init,kernel_regularizer=reg),
    BatchNormalization(),
    Dropout(drop),
    #Dense(8, activation='relu', kernel_initializer=kernel_init, bias_initializer=bias_init,kernel_regularizer=reg),
    #BatchNormalization(),
    #Dropout(drop),
    Dense(1, activation='relu', kernel_initializer=kernel_init, bias_initializer=bias_init,kernel_regularizer=reg)
])
model.compile(optimizer = opt,
              loss = 'mean_squared_error'
              )
history = model.fit(x_train, y_train, epochs = 25, batch_size = 512, validation_split = 0.15, callbacks=[callback])

metrics = pd.DataFrame(history.history)
#metrics.to_excel("modeln.xlsx")
#%%

import matplotlib.pyplot as plt
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Loss vs. epochs')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Training', 'Validation'], loc='upper right')
plt.savefig('modeln.png')
plt.show()
plt.close()
print('Done')

from sklearn.metrics import mean_squared_error
y_train_pred = model.predict(X_train[...,np.newaxis])
print(mean_squared_error(y_train, y_train_pred))

#print(model.evaluate(X_test[...,np.newaxis], y_test))
#%%
# Mice data
mice_data = pd.read_excel('mice_data_MMF.xlsx')
X_mice = mice_data.iloc[:,2:]
X_mice = X_mice.to_numpy()
X_mice = X_mice/100
X_mice = np.reshape(X_mice, (X_mice.shape[0], num_mid, num_metab),'F')
mice_pred = model.predict(X_mice[...,np.newaxis])

pred_results = mice_data.iloc[:,0:2]
pred_results['de_novo_GMP_ratio'] = mice_pred
pred_results.to_excel("MMF_pred.xlsx")

#diff = abs(y_train.to_numpy().reshape(y_train.shape[0],1)-y_train_pred)
#print(sum(diff <= 0.1)/len(y_train_pred))

#%%

def cv_optuna(trial):
    neurons = trial.suggest_int('neurons', low = 8, high = 16, step = 4)
    layers = trial.suggest_int('layers', low = 1, high = 3, step = 1)
    
    learning_rate = trial.suggest_float('learning_rate', low = 0.001, high = 0.01, log = True)
    epsilon = trial.suggest_float('epsilon', low = 0.001, high = 1, log = True)
    #optimizer = Adam(learning_rate = learning_rate, epsilon = epsilon)
    
    #activation = trial.suggest_categorical('activation', ['relu', 'sigmoid', 'tanh'])
    #if activation == 'relu':
    #    kernel_initializer = 'he_uniform'
    #else:
        #kernel_initializer = 'glorot_uniform'
    
    bias = trial.suggest_float('bias', low = 0, high = 0.1, log = True)    
    bias_initializer = tfk.initializers.Constant(value = bias)
    
    reg = trial.suggest_float('reg', low = 0.001, high = 0.1, log = True)
    kernel_regularizer = tfk.regularizers.L2(l2 = reg)
    
    filters = trial.suggest_int('filters', low = 4, high = 16, step = 2)
    
    dropout_rate = trial.suggest_float('dropout_rate', low = 0.1, high = 0.3, step = 0.1)
    
    batch_size = trial.suggest_int('batch_size', low = 64, high = 512, step = 64)
    #epochs = trial.suggest_int('epochs', low = 100, high = 200, step = 50)

    def create_model():
        model = tfk.Sequential()
        model.add(tfk.layers.Conv2D(filters, (no_metabs, 1), 
                                    activation = 'relu', 
                                    kernel_initializer = 'he_uniform', 
                                    kernel_regularizer = kernel_regularizer,
                                    bias_initializer = bias_initializer,
                                    input_shape = (no_metabs, no_mids, 1)))
        model.add(tfk.layers.Flatten())
        model.add(tfk.layers.BatchNormalization())
        model.add(tfk.layers.Dropout(dropout_rate))
        for _ in range(layers):
            model.add(tfk.layers.Dense(neurons, activation = 'relu', 
                                       kernel_initializer = 'he_uniform',
                                       kernel_regularizer = kernel_regularizer,
                                       bias_initializer = bias_initializer))
            model.add(tfk.layers.BatchNormalization())        
            model.add(tfk.layers.Dropout(dropout_rate))
        
        model.add(tfk.layers.Dense(1, kernel_initializer = 'glorot_uniform', 
                                   kernel_regularizer = kernel_regularizer)) 
                                   
        model.compile(loss = 'mean_squared_error', 
                      optimizer = Adam(learning_rate = learning_rate, epsilon = epsilon), 
                      metrics = ['mse'])
        print(model.summary())
        return model
    estimator = KerasRegressor(model = create_model, batch_size = batch_size, 
                               epochs = 200, verbose = 0)
    cv_gfk = GroupKFold(n_splits = 5) # 5-fold split for cross validation   
    es = tfk.callbacks.EarlyStopping(monitor='val_loss', patience = 5, 
                                     start_from_epoch = 50)
                                 #restore_best_weights = True)
    with parallel_backend(backend='threading', n_jobs=-1):
        cv_score = cross_val_score(estimator, train_x, train_y, 
                                   scoring = 'neg_mean_squared_error', 
                                   n_jobs = -1, fit_params={'callbacks':[es]},
                                   groups = train_ind, cv = cv_gfk).mean()
    return cv_score


nn_bo = optuna.create_study(sampler = optuna.samplers.TPESampler(), 
                            direction = 'maximize')
start = time.time()
nn_bo.optimize(cv_optuna, n_trials = 50)
print('It takes %s minutes' % ((time.time() - start)/60))