# Import libraries ------------------------------------------------------------
library(Seurat)
library(tidyr)
library(dplyr)
library(stringr)

# imputation ------------------------------------------------------------------

## Read purine simulated MIDs ----
mid = read.csv(file = paste('../metabolic_CNN/data/sim_data/', 
                            'simulated_data_gmp_denovo_glioma.csv', 
                            sep = ''), header = T, quote = '')


## Read simulation parameters including serine and R5P ----
param_mid = read.csv(file =  paste('../metabolic_CNN/data/sim_data/', 
                                   'simulation_parameters_gmp_denovo_glioma.csv', 
                                   sep = ''), header = T, quote = '')

param_mid[['index']] = seq(1, dim(param_mid)[1])

add_mid <- function(mid_df, added_mid, param_df){
  mid_df[[added_mid]] = param_df[match(mid_df$index, param_df$index), added_mid]
  return(mid_df)
}

mid = add_mid(mid, 'SER1', param_mid) 
mid = add_mid(mid, 'SER2', param_mid)  
mid = add_mid(mid, 'SER3', param_mid)  

mid = add_mid(mid, 'R5P1', param_mid) 
mid = add_mid(mid, 'R5P2', param_mid)  
mid = add_mid(mid, 'R5P3', param_mid)  
mid = add_mid(mid, 'R5P4', param_mid) 
mid = add_mid(mid, 'R5P5', param_mid)  

mid[['SER0']] = 1 - mid$SER1 - mid$SER2 - mid$SER3
mid[['R5P0']] = 1 - mid$R5P1 - mid$R5P2 - mid$R5P3 - mid$R5P4 - mid$R5P5

mid = mid[, !(colnames(mid) %in% c('MTHF0', 'MTHF1', 'GLY0', 'GLY1', 'GLY2',
                                   'SER0', 'R5P0', 'IMP0', 'GMP0', 'GDP0',
                                   'INO0', 'GUO0', 'AMP0'))]

mid$time = as.double(mid$time)

mid = subset(mid, time >= 2, select = colnames(mid))

mid_v = mid[, !(colnames(mid) %in% c('index', 'time'))]

## Read patient MIDs sampled from a truncated normal distribution ----
patient_dir = '../metabolic_CNN/data/patient_data/patient_truncated_norm_samples_purine'
patient_files = list.files(path = patient_dir)
mid_names = read.delim('../metabolic_CNN/data/patient_data/mid_name_patient_purine.txt', 
                       sep = '\t', header = F)

## Combine simulated mids with patient mids ----
for (ps in 1:length(patient_files)){
  
  mid_mc_ps = read.delim(paste(patient_dir, '/', patient_files[ps], sep = ''), 
                         header = F, quote = '')
  colnames(mid_mc_ps) = mid_names$V1
  rownames(mid_mc_ps) = paste(rep(paste('P', ps, sep = '')), 
                              seq(1, dim(mid_mc_ps)[1]), sep = '_')
  
  mid_mc_ps = mid_mc_ps[, !colnames(mid_mc_ps) %in% c('SER0', 'R5P0', 
                                                      'IMP0', 'GMP0', 'GDP0',
                                                      'INO0', 'GUO0', 'AMP0')]
  
  mid_v = rbind(mid_mc_ps, mid_v)
}

mid_v = mid_v * 100
mid_v = data.matrix(mid_v)
mid_v = t(mid_v)

# gbm 1e4
sim_mid = mid_v[, 1501:211500]
patient_mid = mid_v[, 1:1500]

## Save simulated MIDs and patient MIDs as Seurat objects
sim_mid = CreateSeuratObject(counts = sim_mid, project = 'sim_data_purine_v5_1e4')
patient_mid = CreateSeuratObject(counts = patient_mid, project = 'patient_mid_mc')

sim_mid = SCTransform(sim_mid) %>% RunPCA()
sim_mid = RunUMAP(sim_mid, dims = 1:30)

patient_mid = SCTransform(patient_mid) %>% RunPCA()

## Impute patient MIDs into simulated MIDs
anchors = FindTransferAnchors(reference = sim_mid, 
                              query = patient_mid, 
                              normalization.method = "SCT")

predictions.assay = TransferData(anchorset = anchors, 
                                 refdata = GetAssayData(sim_mid[['RNA']]),
                                 prediction.assay = T,
                                 weight.reduction = patient_mid[["pca"]],
                                 dims = 1:30)

pred_assay = as.data.frame(t(as.matrix(predictions.assay@data)))

write.csv(pred_assay, row.names = T, 
          file = 'patient_purine_gbm_2405_1e4_rm0_100.csv')

for (i in 1:15){
  f = (i - 1) * 100 + 1
  e = i * 100
  write.csv(pred_assay[f:e, ], row.names = T,
            file = paste(substr(rownames(pred_assay)[i * 100], 1, 3),
                         'patient_purine_gbm_2405_1e4_rm0_100.csv', sep = '_'))
}

saveRDS(sim_mid, 'sim_glioma_purine_mid_v5_1e4_rm0_100_sct.rds')
saveRDS(patient_mid, 'sim_patient_purine_mid_mc_sct.rds')

## Read MMF-treated mice MIDs sampled from a truncated normal distribution ----
mice_dir = '../metabolic_CNN/data/mice_data/mice_mmf_truncated_norm_samples_purine'
mice_files = list.files(path = mice_dir)
mid_names = read.delim('../metabolic_CNN/data/mice_data/mid_name_mice_purine.txt', 
                       sep = '\t', header = F)

mice_mid = c()
for (ps in 1:length(mice_files)){
  
  mid_mc_ps = read.delim(paste(mice_dir, '/', mice_files[ps], sep = ''), 
                         header = F, quote = '')
  colnames(mid_mc_ps) = mid_names$V1
  rownames(mid_mc_ps) = paste(rep(paste('P', ps, sep = '')), 
                              seq(1, dim(mid_mc_ps)[1]), sep = '_')
  
  mid_mc_ps = mid_mc_ps[, !colnames(mid_mc_ps) %in% c('SER0', 'R5P0', 
                                                      'IMP0', 'GMP0', 'GDP0',
                                                      'INO0', 'GUO0', 'AMP0')]
  
  mice_mid = c(mice_mid, list(mid_mc_ps))
}
names(mice_mid) = substr(mice_files, 1, 7)
mice_mid_v = dplyr::bind_rows(mice_mid)

mice_mid_v = mice_mid_v * 100
mice_mid_v = data.matrix(mice_mid_v)
mice_mid_v = t(mice_mid_v)

mice_mid_v = CreateSeuratObject(counts = mice_mid_v, project = 'mice_mid_mc')

saveRDS(mice_mid_v, 'sim_mice_mid_mc_2405_sct.rds')

mice_mid_v = SCTransform(mice_mid_v) %>% RunPCA()

## Impute mice MIDs into simulated MIDs
anchors = FindTransferAnchors(reference = sim_mid, 
                              query = mice_mid_v, 
                              normalization.method = "SCT")

predictions.assay = TransferData(anchorset = anchors, 
                                 refdata = GetAssayData(sim_mid[['RNA']]),
                                 prediction.assay = T,
                                 weight.reduction = mice_mid_v[["pca"]],
                                 dims = 1:30)

pred_mice_assay = as.data.frame(t(as.matrix(predictions.assay@data)))
write.csv(pred_mice_assay, row.names = T, 
          file = 'mice_purine_gbm_2405_1e4_rm0_100.csv')

for (i in 1:20){
  f = (i - 1) * 100 + 1
  e = i * 100
  write.csv(pred_mice_assay[f:e, ], row.names = T,
            file = paste(names(mice_mid)[i],
                         'mice_purine_gbm_2405_1e4_rm0_100.csv', sep = '_'))
}

## GBM38 TRP

sim_mid = readRDS('./data/sim_data/sim_gbm_purine_mid_v5_1e4_rm0_100_sct.rds')

# mice MIDs
mice_dir = './data/mice_data/mice_mid_mc_purine/'
mice_files = list.files(path = mice_dir)
mid_names = read.delim('./data/mice_data/mid_name_mice_purine_v2.txt', 
                       sep = '\t', header = F)

mice_mid = c()
for (ps in 1:length(mice_files)){
  
  mid_mc_ps = read.delim(paste(mice_dir, 
                               mice_files[ps], sep = ''), 
                         header = F, quote = '')
  colnames(mid_mc_ps) = mid_names$V1
  rownames(mid_mc_ps) = paste(substr(mice_files[ps], 1, nchar(mice_files[ps]) - 27), 
                              seq(1, dim(mid_mc_ps)[1]), sep = '_')
  
  mid_mc_ps = mid_mc_ps[, !colnames(mid_mc_ps) %in% c('ser0', 'r5p0', 
                                                      'IMP0', 'GMP0', 'GDP0',
                                                      'INO0', 'GUO0', 'AMP0')]
  
  mice_mid = c(mice_mid, list(mid_mc_ps))
}
names(mice_mid) = substr(mice_files, 1, nchar(mice_files) - 27)
mice_mid_v = dplyr::bind_rows(mice_mid)

mice_mid_v = mice_mid_v * 100
mice_mid_v = data.matrix(mice_mid_v)
mice_mid_v = t(mice_mid_v)

mice_mid_v = CreateSeuratObject(counts = mice_mid_v, project = 'mice_mid_mc')

saveRDS(mice_mid_v, 'sim_mice_mid_mc_2405_sct.rds')

mice_mid_v = SCTransform(mice_mid_v) %>% RunPCA()
anchors = FindTransferAnchors(reference = sim_mid, 
                              query = mice_mid_v, 
                              normalization.method = "SCT")

predictions.assay = TransferData(anchorset = anchors, 
                                 refdata = GetAssayData(sim_mid[['RNA']]),
                                 prediction.assay = T,
                                 weight.reduction = mice_mid_v[["pca"]],
                                 dims = 1:30)

pred_mice_assay = as.data.frame(t(as.matrix(predictions.assay@data)))

for (i in 1:2){ #20
  f = (i - 1) * 100 + 1
  e = i * 100
  write.csv(pred_mice_assay[f:e, ], row.names = T, quote = F,
            file = paste(names(mice_mid)[i],
                         'purine_gbm_imputed.csv', sep = '_'))
}

