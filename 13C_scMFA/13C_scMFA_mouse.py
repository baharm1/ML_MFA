"""
13C-scMFA - an integrated single cell flux balance analysis and metabolic flux analysis.
The following code integrates scRNA-seq data with MIDs and quantifies 
single cell metabolic fluxes.
"""
import argparse
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
from util import construct_path
from EarlyStopper import EarlyStopper

pd.options.mode.chained_assignment = None  # default='warn'

# hyper parameters
LEARN_RATE = 0.008

def myLoss(
    gene_flux, compound_acc_gene, isotopologue_acc, penalized_mid_bal,
    cellnames_types, geneScale, config_dict):

    # flux balance equations, no accumaltion of compounds in cells
    celltype_counts = cellnames_types['cell_type'].value_counts()
    total1 = torch.pow(compound_acc_gene, 2)
    total1 = torch.sum(total1, dim = 1)
    
    # flux balance equations, no accumaltion of compounds in TME 
    SERsecretion_gene = 0
    # isotopologue balance equations, no accumulation of metabolites in tissue
    total5 = 0
    
    n_mids = config_dict["n_mids"]  
    n_metabs = config_dict["n_balanced_metabs"]
    lamb1 = config_dict["la1_comp_bal"]
    lamb2 = config_dict["la2_non_neg"]
    lamb3 = config_dict["la3_cell_var"]
    lamb4 = config_dict["la4_mid_bal"]
    lamb5 = config_dict["la5_mid_bulk"]
    
    if config_dict['model_name'] == 'serine':
        SERsecretion_gene = (torch.sum(gene_flux[:, 3]) / 
                             celltype_counts['Astrocyte'] - 
                             torch.sum(gene_flux[:, 4]) / 
                             (celltype_counts['Neuron_1'] + 
                              celltype_counts['Neuron_2']) -
                             torch.sum(gene_flux[:, 5]) / 
                             (celltype_counts['Invading_Neoplastic'] +
                             celltype_counts['Neoplastic']))
        SERsecretion_gene = torch.pow(SERsecretion_gene, 2) 
        
        for mid in range(n_mids * n_metabs):
            total5 = total5 + torch.pow(torch.sum(isotopologue_acc[:, mid]) / 
                                        celltype_counts['Astrocyte'] +
                                        torch.sum(isotopologue_acc[:, 
                                                                   mid +  
                                                                   n_mids * 
                                                                   n_metabs]) / 
                                        (celltype_counts['Neuron_1'] + 
                                         celltype_counts['Neuron_2']), 2)
            
            total5 = total5 + torch.pow(torch.sum(isotopologue_acc[:, 
                                                                   mid + 2 *
                                                                   n_mids * 
                                                                   n_metabs]) / 
                                        (celltype_counts['Invading_Neoplastic'] +
                                        celltype_counts['Neoplastic']), 2)
    
    elif config_dict['model_name'] == 'purine':
        for mid in range(n_mids * n_metabs):
            total5 = total5 + torch.pow(torch.sum(isotopologue_acc[:, mid]) / 
                                        celltype_counts['Myeloid'] +
                                        torch.sum(isotopologue_acc[:, mid + 
                                                                   n_mids * 
                                                                   n_metabs]) / 
                                        celltype_counts['Neoplastic'], 2)  
    
    # non-negative constraint
    error = torch.abs(gene_flux) - gene_flux
    total2 = torch.sum(error, dim = 1)
    
    # cell-wise variation constraintt
    total3 = torch.pow(torch.sum(gene_flux, dim = 1) - geneScale, 2) 
    
    # MID accumaltion
    total4 = torch.pow(penalized_mid_bal, 2)
    total4 = torch.sum(total4, dim = 1)

    # loss
    loss1 = torch.sum(lamb1 * total1) + lamb1 * SERsecretion_gene
    loss2 = torch.sum(lamb2 * total2)
    loss3 = torch.sum(lamb3 * total3)
    loss4 = torch.sum(lamb4 * total4)
    loss5 = torch.sum(lamb5 * total5) 
    loss = loss1 + loss2 + loss3 + loss4 + loss5
    
    return loss, loss1, loss2, loss3, loss4, loss5


def main(args: argparse.Namespace):
    #%% set arguments    
    # change working directory to source file location
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    
    # config file
    config_dir = args.config_file
    
    if not isinstance(config_dir, Path):
        config_dir = Path(config_dir)
    with open(config_dir, 'r') as j:
        config_dict = json.loads(j.read())
    
    if not os.path.exists(config_dict['model_dir']):
        raise RuntimeError("Model directory not found.")
            
    data_path = construct_path(config_dict['model_dir'])
    patient_path = construct_path(config_dict['input_data_dir'])
    patient_scRNA_data_files = os.listdir(patient_path / 'scRNA')
    patient_MID_data_files = os.listdir(patient_path / 'MID')
    patient_celltype_files = os.listdir(patient_path / 'celltypes')
    
    moduleGene_file = data_path / config_dict['module_gene_file']
    cm_file = data_path / config_dict['stoichiometry_matrix']
    mm_file = data_path / config_dict['isotopologue_matrix']
    cName_file = config_dict['compound_name_file']
    mName_file = config_dict['isotopologue_name_file']
    
    output_dir = config_dict['output_dir']
    if not os.path.exists(output_dir):
      os.mkdir(output_dir)
      print("Folder %s created!" % output_dir)
    else:
      print("Folder %s already exists" % output_dir)
      
    res_dir = construct_path(config_dict['output_dir'])
    n_mids = config_dict["n_mids"]  
    n_metabs = config_dict["n_balanced_metabs"]
    
    sc_imputation = True
    
    cName = pd.read_csv(data_path / cName_file, sep = ",", header = 0)
    cName = cName.columns
        
    moduleGene = pd.read_csv(moduleGene_file, sep = ",", index_col = 0)
    moduleLen = moduleGene.notnull().sum(axis=1).to_numpy()
    
    # find existing gene
    module_gene_all = moduleGene.to_numpy().flatten()
    module_gene_all = pd.Series(module_gene_all).dropna().to_numpy()
    module_gene_all = set(module_gene_all)
    
    for p in range(len(patient_scRNA_data_files)):
        gene_exp_file = patient_path / 'scRNA' / patient_scRNA_data_files[p]
        cellname_type_file = patient_path / 'celltypes' / patient_celltype_files[p]
        patient_enrich_file = patient_path / 'MID' / patient_MID_data_files[p]
        
        print(patient_MID_data_files[p])
        
        # Sanity checks
        for f in [gene_exp_file, moduleGene_file, cm_file]:
            if not f.is_file():
                raise FileNotFoundError(f"{f} is not found")
        
        # choose cpu or gpu automatically
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        
        # read data
        print("Start loading data ...")
        geneExpr = pd.read_csv(gene_exp_file, index_col = 0)
        geneExpr = geneExpr.T
            
        if sc_imputation == True:
            magic_operator = magic.MAGIC()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                geneExpr = magic_operator.fit_transform(geneExpr)
        if geneExpr.max().max() > 50:
            geneExpr = (geneExpr + 1).apply(np.log2)
        
        if config_dict['model_name'] == "serine":
            ULK_complex = geneExpr[['ATG13_a', 'ATG101_a', 'RB1CC1_a', 
                                    'ATG13_n', 'ATG101_n', 'RB1CC1_n', 
                                    'ATG13_g', 'ATG101_g', 'RB1CC1_g']]
        
            ULK_complex['ULK1_2_a'] = geneExpr['ULK1_a'] + geneExpr['ULK2_a']
            ULK_complex['ULK1_2_n'] = geneExpr['ULK1_n'] + geneExpr['ULK2_n']
            ULK_complex['ULK1_2_g'] = geneExpr['ULK1_g'] + geneExpr['ULK2_g']
        
            geneExpr['ULK_complex_a'] = ULK_complex[['ULK1_2_a','ATG13_a', 
                                                     'ATG101_a', 'RB1CC1_a']].min(axis=1)
            geneExpr['ULK_complex_n'] = ULK_complex[['ULK1_2_n','ATG13_n', 
                                                     'ATG101_n', 'RB1CC1_n']].min(axis=1)
            geneExpr['ULK_complex_g'] = ULK_complex[['ULK1_2_g','ATG13_g', 
                                                     'ATG101_g', 'RB1CC1_g']].min(axis=1)
        
        
        cellnames_types = pd.read_csv(cellname_type_file, sep = '\t', 
                                      index_col = 0, header = 0)
        
        geneExprSum = geneExpr.sum(axis = 1)
        stand = geneExprSum.mean()
        geneExprScale = geneExprSum / stand
        geneExprScale = torch.FloatTensor(geneExprScale.values).to(device)
        
        BATCH_SIZE = geneExpr.shape[0]
        
        data_gene_all = set(geneExpr.columns)
        gene_overlap = list(data_gene_all.intersection(module_gene_all))
        gene_overlap.sort()
        
        cmMat = pd.read_csv(cm_file, sep = ",", header = None)
        cmMat = cmMat.values
        cmMat = torch.FloatTensor(cmMat).to(device)
        
        print("Loading scRNA-seq data is done.")
        
        print("Start processing scRNA-seq data ...")
        emptyNode = []
        # extract overlap gene
        geneExpr = geneExpr[gene_overlap]
        gene_names = geneExpr.columns
        gene_names_set = set(gene_names)
        cell_names = geneExpr.index.astype(str)
        geneExpr = geneExpr.T
        n_modules = moduleGene.shape[0]
        n_genes, n_cells = geneExpr.shape
        
        geneExprDf = pd.DataFrame(columns=["Module_Gene"] + list(cell_names))
        
        for i in range(n_modules):
            genes = set(moduleGene.iloc[i, :].dropna())
            if not genes:
                emptyNode.append(i)
                continue
            temp = geneExpr.copy()
            temp.loc[list(gene_names_set - genes), :] = 0 
            temp["Module_Gene"] = ["%02d_%s" % (i, g) for g in gene_names]
            geneExprDf = pd.concat([geneExprDf, temp], 
                                   ignore_index = True, sort = False)
        geneExprDf.set_index("Module_Gene", inplace = True)
        X = geneExprDf.values.T.astype(float)
        X = torch.FloatTensor(X).to(device)
        
        # prepare data for constraint of module variation based on gene
        df = geneExprDf
        df.index = [i.split("_")[0] for i in df.index]
        
        # change type to ensure correct order
        df.index = df.index.astype(int)
        module_scale = df.groupby(df.index).sum().T
        module_scale = torch.FloatTensor(module_scale.values / moduleLen)
        print("Processing scRNA-seq data is done.")
        
        print("Start loading MID data ...")
        patient_MID = pd.read_csv(patient_enrich_file, sep='\t', index_col=0)
        
        if config_dict['model_name'] == 'serine':
            patient_MID = patient_MID.T
            astro_MIDs = patient_MID[['PGc0', 'PGc1', 'PGc2', 'PGc3', 
                                      'SERc0', 'SERc1', 'SERc2', 'SERc3']].copy()
            astro_MIDs = astro_MIDs.rename(columns={'PGc0': 'PGa0', 'PGc1': 'PGa1', 
                                                    'PGc2': 'PGa2', 'PGc3': 'PGa3',
                                                    'SERc0': 'SERa0', 'SERc1': 'SERa1', 
                                                    'SERc2': 'SERa2', 'SERc3': 'SERa3'})
            neuron_MIDs = patient_MID[['PGc0', 'PGc1', 'PGc2', 'PGc3', 
                                       'SERc0', 'SERc1', 'SERc2', 'SERc3']].copy()
            neuron_MIDs = neuron_MIDs.rename(columns={'PGc0': 'PGn0', 'PGc1': 'PGn1',
                                                      'PGc2': 'PGn2', 'PGc3': 'PGn3',
                                                      'SERc0': 'SERn0', 'SERc1': 'SERn1', 
                                                      'SERc2': 'SERn2', 'SERc3': 'SERn3'})
            
            neoplastic_MIDs = patient_MID[['PGg0', 'PGg1', 'PGg2', 'PGg3', 
                                           'SERg0', 'SERg1', 'SERg2', 'SERg3']].copy()
            plasma_MIDs = patient_MID[['SERp0', 'SERp1', 'SERp2', 'SERp3']].copy()
            
            patient_MID = pd.concat([astro_MIDs, neuron_MIDs, 
                                     neoplastic_MIDs, plasma_MIDs], 
                                    axis = 1)
            patient_MID = patient_MID.T
            
        patient_MID = patient_MID / 100
        
        mmMat = pd.read_csv(mm_file, sep = ",", header = None)
        mmMat = mmMat.astype(str)
        uniq_MID = np.unique(mmMat.to_numpy()).tolist()
        if '0' in uniq_MID: uniq_MID.remove('0')
        if '1' in uniq_MID: uniq_MID.remove('1')
        
        for MID_i in uniq_MID:
            if MID_i[0:3] == 'neg':
                mmMat = mmMat.replace(MID_i, -1 * patient_MID.loc[MID_i[4:], 'mean'])
            else:
                mmMat = mmMat.replace(MID_i, patient_MID.loc[MID_i, 'mean'])
        
        mmMat = mmMat.astype(float)
        mmMat = mmMat.values
        mmMat = torch.FloatTensor(mmMat).to(device)
        
        mName = pd.read_csv(data_path / mName_file, sep = ",", header = 0)
        mName = mName.columns
        
        print("Loading MID data is done.")
            
        #%% Initialize a neural network and optimize it by minimizing obj func
        torch.manual_seed(16)
        net = FLUX(n_modules, n_genes).to(device)
        
        optimizer = torch.optim.Adam(net.parameters(), lr=LEARN_RATE)
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
        train_loader = DataLoader(dataset=dataSet, **dataloader_params)
        
        # Train neural network
        print("Start training neural network ...")
        start = time.time()
        # save loss values for 5 terms of obj func
        loss_v = []
        loss_v1 = []
        loss_v2 = []
        loss_v3 = []
        loss_v4 = []
        loss_v5 = []
        
        net.train()
        timestr = time.strftime("%Y%m%d-%H%M%S")
        lossName = res_dir / f"lossValue_{timestr}.txt"
        file_loss = open(lossName, "a")
        epoch = 0
        
        def generator():
          while True:
            yield
        
        
        for _ in tqdm(generator()):
        
            epoch += 1
        
            
            loss, loss1, loss2, loss3, loss4, loss5 = 0, 0, 0, 0, 0, 0
        
            for i, (X, X_scale, m_scale) in enumerate(train_loader):
        
                X_batch = Variable(X.float().to(device))
                X_scale_batch = Variable(X_scale.float().to(device))
                m_scale_batch = Variable(m_scale.float().to(device))
                
                
                out_mg_batch, out_cg_batch, out_ig_batch, out_mid_batch = net(
                    X_batch, m_scale_batch, cmMat, mmMat)
                
                
                (loss_batch, loss1_batch, loss2_batch, loss3_batch, loss4_batch, 
                 loss5_batch) = myLoss(
                    out_mg_batch, out_cg_batch, out_ig_batch, out_mid_batch,
                    cellnames_types, 
                    geneScale = X_scale_batch,
                    config_dict = config_dict)
        
                optimizer.zero_grad() # reset the gradients to 0
                loss_batch.backward() # perform a backward pass
                optimizer.step() # update the weights
         
                loss += loss_batch.cpu().data.numpy()
                loss1 += loss1_batch.cpu().data.numpy()
                loss2 += loss2_batch.cpu().data.numpy()
                loss3 += loss3_batch.cpu().data.numpy()
                loss4 += loss4_batch.cpu().data.numpy()
                loss5 += loss5_batch.cpu().data.numpy()
                
            
            loss_v.append(loss)
            loss_v1.append(loss1)
            loss_v2.append(loss2)
            loss_v3.append(loss3)
            loss_v4.append(loss4)
            loss_v5.append(loss5)
            
            scheduler.step(loss)
            
            file_loss.write(
                "epoch: %02d, loss1: %.8f, loss2: %.8f, loss3: %.8f, loss4: %.8f, , loss5: %.8f, loss: %.8f, lr: %.8f. \n"
                % (epoch, loss1, loss2, loss3, loss4, loss5, loss, scheduler.optimizer.param_groups[0]['lr'])
            )
            print("epoch: %02d, loss1: %.4f, loss2: %.4f, loss3: %.4f, loss4: %.4f, , loss5: %.4f, loss: %.4f, lr: %.4f."
            % (epoch, loss1, loss2, loss3, loss4, loss5, loss, scheduler.optimizer.param_groups[0]['lr']))
            
            if early_stopper.early_stop(loss):             
                break 
                 
        end = time.time()
        print("Training time (minutes): ", (end - start) / 60)
        
        file_loss.close()
        
        # Save output of neural netwrok: fluxes, model hyperparameters, and balances
        plt.plot(loss_v, "--", c='#4053D3')
        plt.plot(loss_v1, c='#DDB310')
        plt.plot(loss_v2, c='#B51D14')
        plt.plot(loss_v3, c='#CACACA')
        plt.plot(loss_v4, c='#00BEFF')
        plt.plot(loss_v5, c='#FB49B0')
        
        plt.legend(["total", "L1", "L2", "L3", "L4", "L5"])
        imgName = res_dir / f"loss_{timestr}.pdf"
        plt.savefig(imgName, format = 'pdf')
        plt.close()
        
        model_params = pd.DataFrame(
            {'run_time': end - start, 'compoundBalance_g': config_dict["la1_comp_bal"], 
             'negativeFlux_g': config_dict["la2_non_neg"], 'cellVar_g': config_dict["la3_cell_var"],
             'moduleVar_g': config_dict["la4_mid_bal"], 'isotoplogueBalance_g': config_dict["la5_mid_bulk"], 
             'n_modules': n_modules, 'n_metabs': n_metabs, 'n_mids': n_mids, 
             'n_cells': n_cells, 'n_batches': BATCH_SIZE, 'n_epochs': epoch,
             'sc_imputation': sc_imputation}, index = [0]).T
        
        model_params.to_csv(res_dir / f"model_params_{timestr}.txt", sep = '\t')
        
        # Run trained model on data again
        net.eval()
        
        for i, (X, X_scale, m_scale) in enumerate(train_loader):
        
            X_batch = Variable(X.float().to(device))
            m_scale_batch = Variable(m_scale.float().to(device))
            
            out_mg_batch, out_cg_batch, out_ig_batch, out_mid_batch = net(
                X_batch, m_scale_batch, cmMat, mmMat)
        
            fluxBalance_g = out_mg_batch.detach().cpu().numpy() 
            compBalance_g = out_cg_batch.detach().cpu().numpy()
            midBalance_g = out_mid_batch.detach().cpu().numpy()
        
        # save to file
        setFg = pd.DataFrame(fluxBalance_g)
        setFg.columns = moduleGene.index
        setFg.index = geneExpr.columns
        setFg.to_csv(res_dir / (f"{gene_exp_file.stem}_fg_{timestr}.csv"))
        
        setBg = pd.DataFrame(compBalance_g)
        setBg.rename(columns=lambda x: x + 1)
        setBg.index = geneExpr.columns
        if cName_file != "noCompoundName":
            setBg.columns = cName
        setBg.to_csv(res_dir / f"balance_cg_{timestr}.csv")
        
        
        setIm = pd.DataFrame(midBalance_g)
        setIm.rename(columns=lambda x: x + 1)
        setIm.index = geneExpr.columns
        if mName_file != "noCompoundName":
            setIm.columns = mName
        setIm.to_csv(res_dir / f"balance_im_{timestr}.csv")
        
        print("13C-scMFA job finished successfully.")
        
    return


def parse_arguments(parser: argparse.ArgumentParser):

    parser.add_argument(
        "--config_file",
        type = str,
        default = "./serine_13C_scMFA_mouse_config.json",
        help = "Input arguments of the model are provided in a JSON file which contains the following information: \n"
        "model_name: Name of the metabolic model (serine or purine) \n"
	    "model_dir: Directory of the metabolic model \n"
	    "input_data_dir: Directory of mouse data, including mouse MIDs (percent enrichment), scRNA-seq data in the form of a block diagonal gene expression matrix constructed in prep_input_data_scMFA_mouse.R, and cell types assigned to cell IDs in scRNA-seq data. These subdirectories are as follows: MID, scRNA, and celltypes. \n"
	    "output_dir: Directory of <sup>13</sup>C-scMFA output files \n"
	    "module_gene_file: Filename of genes associated with metabolic reactions in model_dir \n"
	    "stoichiometry_matrix: Filename of metabolite balances in model_dir \n"
	    "isotopologue_matrix: Filename of isotopologue balances in model_dir \n"
	    "compound_name_file: Filename of balanced metabolite names in model_dir \n"
	    "isotopologue_name_file: Filename of balanced isotopologue names in model_dir \n"
	    "n_mids: Number of balanced isotopologues for a balanced metabolite \n"
	    "n_balanced_metabs: Number of balanced metabolites in a cell \n"
        "la1_comp_bal: Hyperparameter to adjust accumulation of metabolites \n"
        "la2_non_neg: Hyperparameter to adjust negative fluxes \n" 
        "la3_cell_var: Hyperparameter to adjust single cell flux variation \n"
        "la4_mid_bal: Hyperparameter to adjust accumulation of MIDs in cells \n"
        "la5_mid_bulk: Hyperparameter to adjust accumulation of MIDs in tissues \n"
        "Please see the default JSON file as an example."
    )
    
    arg = parser.parse_args()

    return arg

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="13C-scMFA - an integrated single cell flux balance analysis and metabolic flux analysis",
        formatter_class=argparse.RawTextHelpFormatter
    )
    args = parse_arguments(parser)
    main(args)    