"""
The following code is a modified version of scFEA 
(https://github.com/changwn/scFEAAlghamdi et al, Genome Research, 2021) 
to include exchange reactions between cell types.
"""
import json
import os
from pathlib import Path
import time
import warnings
import magic
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from torch.autograd import Variable
from torch.utils.data import DataLoader
from tqdm import tqdm

from ClassFlux import FLUX  # Flux class network
from DatasetFlux import MyDataset
from util import pearsonr
from EarlyStopper import EarlyStopper

# hyper parameters
LEARN_RATE = 0.008

LAMB_BA = 1.5
LAMB_NG = 1
LAMB_CELL = 1
LAMB_MOD = 1e-2

#%% Objective function
def myLoss(
    m, c, cellnames_types, lamb1=0.2, lamb2=0.2, lamb3=0.2, lamb4=0.2, 
    geneScale=None, moduleScale=None):
    
    # flux balance equations, no accumulation of compounds
    celltype_counts = cellnames_types['cell_type'].value_counts()
    total1 = torch.pow(c, 2)
    total1 = torch.sum(total1, dim = 1)
    SERsecretion = (torch.sum(m[:, 3]) / celltype_counts['Astrocyte'] - 
                    torch.sum(m[:, 4]) / celltype_counts['Neoplastic'] - 
                    torch.sum(m[:, 5]) / celltype_counts['Neuron'])
    SERsecretion = torch.pow(SERsecretion, 2)

    # non-negative constrain
    error = torch.abs(m) - m
    total2 = torch.sum(error, dim = 1)
    
    # sample-wise variation constrain
    total3 = torch.pow(torch.sum(m, dim = 1) - geneScale, 2) 
    
    # module-wise variation constrain
    if lamb4 > 0:
        corr = torch.FloatTensor(np.ones(m.shape[0]))
        for i in range(m.shape[0]):
            corr[i] = pearsonr(m[i, :], moduleScale[i, :])
        corr = torch.abs(corr)
        penal_m_var = torch.FloatTensor(np.ones(m.shape[0])) - corr
        total4 = penal_m_var
    else:
        total4 = torch.FloatTensor(np.zeros(m.shape[0]))

    # loss
    loss1 = torch.sum(lamb1 * total1) + lamb1 * SERsecretion
    loss2 = torch.sum(lamb2 * total2)
    loss3 = torch.sum(lamb3 * total3)
    loss4 = torch.sum(lamb4 * total4)
    loss = loss1 + loss2 + loss3 + loss4
    
    return loss, loss1, loss2, loss3, loss4

#%% Set arguments

# change working directory to source file location
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# config file
config_dir = './modified_scFEA_config.json'

if not isinstance(config_dir, Path):
    config_dir = Path(config_dir)
with open(config_dir, 'r') as j:
    config_dict = json.loads(j.read())

if not os.path.exists(config_dict['data_dir']):
    raise RuntimeError("Data directory not found.")
    
data_dir = config_dict['data_dir']
res_dir = config_dict['res_dir']
moduleGene_file = "".join([data_dir, '/', config_dict['module_gene_file']])
cm_file = "".join([data_dir, '/', config_dict['stoichiometry_matrix']])
gene_exp_file = "".join([data_dir, '/', config_dict['gene_exp']])
cellnames_types_file = "".join([data_dir, '/', config_dict['cellname_types']])
sc_imputation = True
cName_file = "".join([data_dir, '/', config_dict['compound_name_file']])

# choose cpu or gpu automatically
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# read data
print("Starting load data...")
geneExpr = pd.read_csv(gene_exp_file, index_col = 0)
geneExpr = geneExpr.T

if sc_imputation == True:
    magic_operator = magic.MAGIC()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        geneExpr = magic_operator.fit_transform(geneExpr)
if geneExpr.max().max() > 50:
    geneExpr = (geneExpr + 1).apply(np.log2)
geneExprSum = geneExpr.sum(axis = 1)
stand = geneExprSum.mean()
geneExprScale = geneExprSum / stand
geneExprScale = torch.FloatTensor(geneExprScale.values).to(device)

BATCH_SIZE = geneExpr.shape[0]

moduleGene = pd.read_csv(moduleGene_file, sep = ",", index_col = 0)
moduleLen = moduleGene.notnull().sum(axis = 1).to_numpy()

# find existing gene
module_gene_all = moduleGene.to_numpy().flatten()
module_gene_all = pd.Series(module_gene_all).dropna().to_numpy()
module_gene_all = set(module_gene_all)

data_gene_all = set(geneExpr.columns)
gene_overlap = list(data_gene_all.intersection(module_gene_all))  # fix
gene_overlap.sort()

cellnames_types = pd.read_csv(cellnames_types_file, sep = '\t', index_col = 0,
                              header = 0)

cmMat = pd.read_csv(cm_file, sep = ",", header = None)
cmMat = cmMat.values
cmMat = torch.FloatTensor(cmMat).to(device)

cName = pd.read_csv(cName_file, sep = ",", header = 0)
cName = cName.columns
print("Load data done.")

print("Starting process data...")
emptyNode = []
# extract overlap gene
geneExpr = geneExpr[gene_overlap]
gene_names = geneExpr.columns
gene_names_set = set(gene_names)
cell_names = geneExpr.index.astype(str)
geneExpr = geneExpr.T
n_modules = moduleGene.shape[0]
n_genes, n_cells = geneExpr.shape
n_comps = cmMat.shape[0]
geneExprDf = pd.DataFrame(columns=["Module_Gene"] + list(cell_names))
for i in range(n_modules):
    genes = set(moduleGene.iloc[i, :].dropna())
    if not genes:
        emptyNode.append(i)
        continue
    temp = geneExpr.copy()
    temp.loc[list(gene_names_set - genes), :] = 0.0
    temp["Module_Gene"] = ["%02d_%s" % (i, g) for g in gene_names]
    geneExprDf = pd.concat([geneExprDf, temp], ignore_index=True, sort=False)
geneExprDf.set_index("Module_Gene", inplace=True)
X = geneExprDf.values.T.astype(float)
X = torch.FloatTensor(X).to(device)

# prepare data for constraint of module variation based on gene
df = geneExprDf
df.index = [i.split("_")[0] for i in df.index]

# change type to ensure correct order
df.index = df.index.astype(int)
module_scale = df.groupby(df.index).sum().T
module_scale = torch.FloatTensor(module_scale.values / moduleLen)
print("Process data done.")

#%% Initialize a neural network and optimize it by minimizing obj func
torch.manual_seed(16)
net = FLUX(n_modules, f_in = n_genes, f_out = 1).to(device)
optimizer = torch.optim.Adam(net.parameters(), lr = LEARN_RATE)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 
                                                       mode = 'min', 
                                                       factor = 0.9, 
                                                       patience = 10, 
                                                       threshold = 0.001,
                                                       min_lr = 0.001,
                                                       cooldown = 10)

early_stopper = EarlyStopper(patience = 100, min_delta = 1e-4)

# Prepare data for neural network
dataloader_params = {
    "batch_size": BATCH_SIZE,
    "shuffle": False,
    "num_workers": 0,
    "pin_memory": False,
}

dataSet = MyDataset(X, geneExprScale, module_scale)
train_loader = DataLoader(dataset = dataSet, **dataloader_params)

# Train neural network
print("Start training neural network ...")
start = time.time()

loss_v = []
loss_v1 = []
loss_v2 = []
loss_v3 = []
loss_v4 = []

net.train()
timestr = time.strftime("%Y%m%d-%H%M%S")
lossName = "".join([res_dir, '/', f"lossValue_{timestr}.txt"])
file_loss = open(lossName, "a")

epoch = 0

def generator():
  while True:
    yield


for _ in tqdm(generator()):

    epoch += 1

    loss, loss1, loss2, loss3, loss4 = 0, 0, 0, 0, 0

    for i, (X, X_scale, m_scale) in enumerate(train_loader):

        X_batch = Variable(X.float().to(device))
        X_scale_batch = Variable(X_scale.float().to(device))
        m_scale_batch = Variable(m_scale.float().to(device))

        out_m_batch, out_c_batch = net(X_batch, n_modules, n_genes, n_comps, cmMat)
        loss_batch, loss1_batch, loss2_batch, loss3_batch, loss4_batch = myLoss(
            out_m_batch,
            out_c_batch,
            cellnames_types = cellnames_types,
            lamb1 = LAMB_BA,
            lamb2 = LAMB_NG,
            lamb3 = LAMB_CELL,
            lamb4 = LAMB_MOD,
            geneScale = X_scale_batch,
            moduleScale = m_scale_batch)

        optimizer.zero_grad() # reset the gradients to 0
        loss_batch.backward() # perform a backward pass
        optimizer.step() # update the weights
 
        loss += loss_batch.cpu().data.numpy()
        loss1 += loss1_batch.cpu().data.numpy()
        loss2 += loss2_batch.cpu().data.numpy()
        loss3 += loss3_batch.cpu().data.numpy()
        loss4 += loss4_batch.cpu().data.numpy()

    loss_v.append(loss)
    loss_v1.append(loss1)
    loss_v2.append(loss2)
    loss_v3.append(loss3)
    loss_v4.append(loss4)
    
    scheduler.step(loss)
    
    file_loss.write(
        "epoch: %02d, loss1: %.8f, loss2: %.8f, loss3: %.8f, loss4: %.8f, , loss: %.8f, lr: %.8f. \n"
        % (epoch, loss1, loss2, loss3, loss4, loss, scheduler.optimizer.param_groups[0]['lr'])
    )
    print("epoch: %02d, loss1: %.4f, loss2: %.4f, loss3: %.4f, loss4: %.4f, , loss: %.4f, lr: %.4f."
    % (epoch, loss1, loss2, loss3, loss4, loss, scheduler.optimizer.param_groups[0]['lr']))
    
    if early_stopper.early_stop(loss):             
        break 

end = time.time()
print("Training time: ", end - start)

file_loss.close()

#%% Save output of neural netwrok: fluxes, model hyperparameters, and balances
plt.plot(loss_v, "--")
plt.plot(loss_v1)
plt.plot(loss_v2)
plt.plot(loss_v3)
plt.plot(loss_v4)
plt.legend(["total", "balance", "negative", "cellVar", "moduleVar"])
imgName = "".join([res_dir, '/', f"loss_{timestr}.pdf"])
plt.savefig(imgName, format = 'pdf')
plt.close()

fluxStatuTest = np.zeros((n_cells, n_modules), dtype = "f")  # float32
balanceStatus = np.zeros((n_cells, n_comps), dtype = "f")

model_params = pd.DataFrame(
    {'run_time': end - start, 'compoundBalance_g': LAMB_BA, 
     'negativeFlux_g': LAMB_NG, 'cellVar_g': LAMB_CELL,
     'moduleVar_g': LAMB_MOD, 
     'n_modules': n_modules, 
     'n_cells': n_cells, 'n_batches': BATCH_SIZE, 'n_epochs': epoch,
     'sc_imputation': sc_imputation}, index = [0]).T

model_params.to_csv("".join([res_dir, '/', f"model_params_{timestr}.txt"]), 
                    sep = '\t')

# Run trained model on data again
net.eval()

for i, (X, X_scale, _) in enumerate(train_loader):

    X_batch = Variable(X.float().to(device))
    out_m_batch, out_c_batch = net(X_batch, n_modules, n_genes, n_comps, cmMat)

    fluxStatuTest = out_m_batch.detach().cpu().numpy()
    balanceStatus = out_c_batch.detach().cpu().numpy() 

# save to file
fileName = "".join([res_dir, '/', f"three_comp_scRNA_f_{timestr}.csv"])
    
setF = pd.DataFrame(fluxStatuTest)
setF.columns = moduleGene.index
setF.index = geneExpr.columns
setF.to_csv(fileName)

setB = pd.DataFrame(balanceStatus)
setB.rename(columns=lambda x: x + 1)
setB.index = setF.index
setB.columns = cName

balanceName = "".join([res_dir, '/', f"balance_{timestr}.csv"])
setB.to_csv(balanceName)

print("Modified scFEA job finished. Check result in output folder.")
