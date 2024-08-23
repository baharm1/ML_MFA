# Import libraries ------------------------------------------------------------
library(tidyr)
library(RColorBrewer)
library(see)
library(ggplot2)
library(dplyr)
library(gghalves)
library(stringr)
library(Seurat)

# Plots for serine simulated MIDs ---------------------------------------------

## Read serine simulated MIDs ----
rm(list = ls())

mid_cortex = read.csv(file = paste('../metabolic_CNN/data/sim_data/', 
                                   'simulated_data_serine_denovo_cortex.csv', 
                                   sep = ''), header = T, quote = '')
mid_denovo = read.csv(file = paste('../metabolic_CNN/data/sim_data/', 
                                   'simulated_data_serine_denovo_glioma.csv',
                                   sep = ''), header = T, quote = '')
mid_plasma = read.csv(file = paste('../metabolic_CNN/data/sim_data/', 
                                   'simulated_data_serine_plasma_glioma.csv',
                                   sep = ''), header = T, quote = '')

mid = rbind(mid_denovo, mid_plasma, mid_cortex)

## Read patient MIDs ----
patient_mids = read.delim(file = 'patient_MIDs_cortex_glioma_serine.txt', 
                          header = T, sep = '\t', quote = '')

patient_mids$INFUSION_TIME = patient_mids$INFUSION_TIME/60 # hour

max(patient_mids$INFUSION_TIME) # 3.66 hr
min(patient_mids$INFUSION_TIME) # 2.2 hr

patient_time_points = unique(round(patient_mids$INFUSION_TIME, digits = 1))
patient_time_points

patient_mid_rm = subset(patient_mids, select = -c(PG0, PGc0, SER0, SERc0, 
                                                  GLY0, GLYc0))

head(patient_mid_rm)

patient_mid_rm$PATIENT = str_sub(patient_mid_rm$PATIENT, start = -1)
patient_mid_rm$SITENO = str_sub(patient_mid_rm$SITE, start = -1)
patient_mid_rm$SITE = ifelse(str_sub(patient_mid_rm$SITE, start = 1, 
                                     end = -2) == 'E', 1, 0)

patient_mid_rm$INFUSION_TIME = round(patient_mid_rm$INFUSION_TIME, digits = 1)

patient_mid_long = pivot_longer(patient_mid_rm, cols = 5:ncol(patient_mid_rm)-1)
head(patient_mid_long)

patient_mid_long[['PATIENTSITE']] = paste(patient_mid_long$PATIENT, 
                                          patient_mid_long$SITE, sep = '')
patient_mid_long[patient_mid_long$PATIENTSITE == 11, 'hexcode'] = '#0e9b74'
patient_mid_long[patient_mid_long$PATIENTSITE == 10, 'hexcode'] = '#6b8740'
patient_mid_long[patient_mid_long$PATIENTSITE == 41, 'hexcode'] = '#447390'
patient_mid_long[patient_mid_long$PATIENTSITE == 40, 'hexcode'] = '#74c4ea'
patient_mid_long[patient_mid_long$PATIENTSITE == 61, 'hexcode'] = '#0087ae'
patient_mid_long[patient_mid_long$PATIENTSITE == 60, 'hexcode'] = '#3fc4e4'
patient_mid_long[patient_mid_long$PATIENTSITE == 21, 'hexcode'] = '#0087ae'
patient_mid_long[patient_mid_long$PATIENTSITE == 20, 'hexcode'] = '#3fc4e4'
patient_mid_long[patient_mid_long$PATIENTSITE == 31, 'hexcode'] = '#0087ae'
patient_mid_long[patient_mid_long$PATIENTSITE == 30, 'hexcode'] = '#3fc4e4'
patient_mid_long[patient_mid_long$PATIENTSITE == 51, 'hexcode'] = '#0087ae'
patient_mid_long[patient_mid_long$PATIENTSITE == 50, 'hexcode'] = '#3fc4e4'
patient_mid_long[patient_mid_long$PATIENTSITE == 71, 'hexcode'] = '#0087ae'
patient_mid_long[patient_mid_long$PATIENTSITE == 70, 'hexcode'] = '#3fc4e4'
patient_mid_long[patient_mid_long$PATIENTSITE == 80, 'hexcode'] = '#0087ae'


sim_mid = mid[mid$time %in% patient_time_points, ]
head(sim_mid)

sim_mid[, seq(3, ncol(sim_mid))] = sim_mid[, seq(3, ncol(sim_mid))] * 100

sim_mid_rmM0 = subset(sim_mid, select = -c(GLY0, GLYc0, SER0, SERc0, 
                                           MTHF0, MTHF1, MTHFc0, MTHFc1, 
                                           PG0, PGc0))
head(sim_mid_rmM0)
ncol(sim_mid_rmM0)

sim_mid_long = pivot_longer(sim_mid_rmM0, cols = 3:ncol(sim_mid_rmM0))
sim_mid_long[['Mi']] = str_sub(sim_mid_long$name, start = -1)
colors = c("#FFC300", "#DA8A67", "#E41B17", "#7D0552", "#2F0909")

for (i in seq(1, 5)){
  sim_mid_long[sim_mid_long$Mi == i, 'colors'] = colors[i]  
}


gc()

## Visualize distribution of MIDs at patient time points ----
for (i in patient_time_points){
  sim_mid_t = sim_mid_long[sim_mid_long$time == i, ]
  patient_points = patient_mid_long[patient_mid_long$INFUSION_TIME == i, ]
  p = ggplot() + 
    
    geom_violinhalf(data = sim_mid_t, aes(x = name, y = value, fill = Mi), 
                    trim = T, lwd = 0.3, scale = 'width') +
    
    geom_half_boxplot(data = sim_mid_t, aes(x = name, y = value),
                      width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
    geom_point(data = patient_points, aes(x = name, y = value, 
                                          colour = PATIENTSITE),
               shape = 1, size = 0.1) + 
    
    theme_classic()+
    labs(
      x    = "Metabolite MID",
      y    = "Simulated MIDs",
      fill = "M+i"
    ) +
    scale_fill_manual(values = colors) +
    coord_flip() + 
    theme(legend.key.size = unit(0.05, 'in'), 
          legend.title = element_text(size = 7))
  
  ggsave(plot = p, filename = paste('sim_MID', i, '.pdf', sep = '_'),
         device = "pdf", height = 4,
         width = 4, units = "in")
  
}

## Visualize distribution of MIDs all timepoints ----
p = ggplot() + 
  geom_violinhalf(data = sim_mid_long, aes(x = name, y = value, fill = Mi), 
                  trim = T, lwd = 0.3, scale = 'width') +
  
  geom_half_boxplot(data = sim_mid_long, aes(x = name, y = value),
                    width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_point(data = patient_mid_long, aes(x = name, y = value, colour = PATIENTSITE),
             shape = 1, size = 0.1) + 
  
  theme_classic()+
  labs(
    x    = "Metabolite MID",
    y    = "Simulated MIDs",
    fill = "M+i"
  ) +
  scale_fill_manual(values = colors) +
  coord_flip() + 
  theme(legend.key.size = unit(0.05, 'in'), legend.title = element_text(size=7))

ggsave(plot = p, filename =  'sim_MID_all_patients_trim_all_glioma_denovo_plasma_cortex_50.pdf',
       device = "pdf", height = 4,
       width = 4, units = "in")

# tsne plots for serine MID simulation ----------------------------------------

## Read simulation parameters including plasma serine ----
param_cortex = read.csv(file = paste('../metabolic_CNN/data/sim_data/', 
                                     'simulation_parameters_serine_denovo_cortex.csv', 
                                     sep = ''), header = T, quote = '')
param_denovo = read.csv(file = paste('../metabolic_CNN/data/sim_data/', 
                                     'simulation_parameters_serine_denovo_glioma.csv',
                                     sep = ''), header = T, quote = '')
param_plasma = read.csv(file = paste('../metabolic_CNN/data/sim_data/', 
                                     'simulation_parameters_serine_plasma_glioma.csv',
                                     sep = ''), header = T, quote = '')

param_plasma[['index']] = seq(1, dim(param_plasma)[1])
param_denovo[['index']] = seq(1, dim(param_denovo)[1])
param_cortex[['index']] = seq(1, dim(param_cortex)[1])

add_mid <- function(mid_df, added_mid, param_df){
  mid_df[[added_mid]] = param_df[match(mid_df$index, param_df$index), added_mid]
  return(mid_df)
}

## Combine plasma serine MIDs from parameters to mid dataframe ----
mid_plasma = add_mid(mid_plasma, 'SERp1', param_plasma) 
mid_plasma = add_mid(mid_plasma, 'SERp2', param_plasma)  
mid_plasma = add_mid(mid_plasma, 'SERp3', param_plasma)  

mid_denovo = add_mid(mid_denovo, 'SERp1', param_denovo) 
mid_denovo = add_mid(mid_denovo, 'SERp2', param_denovo)  
mid_denovo = add_mid(mid_denovo, 'SERp3', param_denovo)  

mid_cortex = add_mid(mid_cortex, 'SERp1', param_cortex) 
mid_cortex = add_mid(mid_cortex, 'SERp2', param_cortex)  
mid_cortex = add_mid(mid_cortex, 'SERp3', param_cortex)  

mid_plasma[['SERp0']] = 1 - mid_plasma$SERp1 - mid_plasma$SERp2 - mid_plasma$SERp3
mid_denovo[['SERp0']] = 1 - mid_denovo$SERp1 - mid_denovo$SERp2 - mid_denovo$SERp3
mid_cortex[['SERp0']] = 1 - mid_cortex$SERp1 - mid_cortex$SERp2 - mid_cortex$SERp3

mid_plasma = mid_plasma[, !(colnames(mid_plasma) %in% c('MTHF0', 'MTHF1', 'MTHFc0', 'MTHFc1'))]
mid_denovo = mid_denovo[, !(colnames(mid_denovo) %in% c('MTHF0', 'MTHF1', 'MTHFc0', 'MTHFc1'))]
mid_cortex = mid_cortex[, !(colnames(mid_cortex) %in% c('MTHF0', 'MTHF1', 'MTHFc0', 'MTHFc1'))]


mid_plasma$time = as.double(mid_plasma$time)

mid_plasma = subset(mid_plasma, time >= 2, select = colnames(mid_plasma))
mid_denovo = subset(mid_denovo, time >= 2, select = colnames(mid_denovo))
mid_cortex = subset(mid_cortex, time >= 2, select = colnames(mid_cortex))

mid_plasma_v = mid_plasma[, !(colnames(mid_plasma) %in% c('index', 'time'))]
mid_denovo_v = mid_denovo[, !(colnames(mid_denovo) %in% c('index', 'time'))]
mid_cortex_v = mid_cortex[, !(colnames(mid_cortex) %in% c('index', 'time'))]

## Read patient MIDs sampled from a truncated normal distribution ----

patient_dir = '../metabolic_CNN/data/patient_data/patient_truncated_norm_samples_serine'
patient_files = list.files(path = patient_dir)

mid_names = read.delim('../metabolic_CNN/data/patient_data/mid_name_patient_serine.txt', 
                       sep = '\t', header = F)

## Combine simulated mids with patient mids ----
for (ps in 1:length(patient_files)){
  
  mid_mc_ps = read.delim(paste(patient_dir, '/', patient_files[ps], sep = ''),
                         header = F, quote = '')
  colnames(mid_mc_ps) = mid_names$V1
  rownames(mid_mc_ps) = paste(rep(paste('P', ps, sep = '')), 
                              seq(1, dim(mid_mc_ps)[1]), sep = '_')
  
  mid_plasma_v = rbind(mid_mc_ps, mid_plasma_v)
  mid_denovo_v = rbind(mid_mc_ps, mid_denovo_v)
  mid_cortex_v = rbind(mid_mc_ps, mid_cortex_v)
}

mid_plasma_v = data.matrix(mid_plasma_v)
mid_plasma_v = t(mid_plasma_v)

com_plasma = CreateSeuratObject(counts = mid_plasma_v, 
                                project = "combined_sim_patient_mid")

com_plasma = FindVariableFeatures(com_plasma, selection.method = "vst")

com_plasma = ScaleData(com_plasma, features = rownames(com_plasma))

com_plasma = RunPCA(com_plasma, 
                    features = VariableFeatures(object = com_plasma), npcs = 10)

ElbowPlot(com_plasma)

com_plasma = FindNeighbors(com_plasma, dims = 1:5)

com_plasma = FindClusters(com_plasma, resolution = 0.1)

com_plasma = RunTSNE(com_plasma, dims = 1:5)

my_colors = c('P8N' = "#5A5156",
              'P7N' = "#C4451C",
              "P7E" = "#822E1C",
              "P6N" = "#1C7F93",
              "P6E" = "#3B00FB",
              "P5N" = "#D85FF7",
              "P5E" = "#FEAF16",
              "P4N" = "#90AD1C",
              "P4E" = "#1C8356",
              "P3N" = "#AA0DFE",
              "P3E" = "#DEA0FD",
              "P2N" = "#2ED9FF",
              "P2E" = "#3283FE",
              "P1N" = "#B00068",
              "P1E" = "#F6222E",
              'sim' = "#E4E1E3")

PS = c(rep('P8N', 100), rep('P7N', 100), rep('P7E', 100), 
       rep('P6N', 100), rep('P6E', 100), rep('P5N', 100), rep('P5E', 100),
       rep('P4N', 100), rep('P4E', 100), rep('P3N', 100), rep('P3E', 100),
       rep('P2N', 100), rep('P2E', 100), rep('P1N', 100), rep('P1E', 100),
       rep('sim', 1050000))

com_plasma@meta.data[['ps_sim']] = PS

pdf(file = 'dimplot_tsne_plasma_patient_site.pdf', height = 10, width = 10)
DimPlot(com_plasma, reduction = "tsne", group.by = 'ps_sim', 
        cells = colnames(com_plasma)[0:1500], cols = my_colors)
dev.off()

pdf(file = 'dimplot_tsne_plasma_sim.pdf', height = 10, width = 10)
DimPlot(com_plasma, reduction = "tsne", group.by = 'ps_sim', 
        cells = colnames(com_plasma)[1501:1051500], cols = my_colors)
dev.off()

saveRDS(com_plasma, 'sim_glioma_50_plasma_ps.rds')

# serine flux distribution ----------------------------------------------------

# Read simulated fluxes
flux_plasma = read.delim(file = './sim_flux/flux_glioma_plasma_serine.txt',
                         header = F, sep = '\t', quote = '')

flux_denovo = read.delim(file = './sim_flux/flux_glioma_denovo_serine.txt',
                         header = F, sep = '\t', quote = '')

flux_cortex = read.delim(file = './sim_flux/flux_cortex_denovo_serine.txt',
                         header = F, sep = '\t', quote = '')

rxn_ids = read.delim(file = './sim_flux/rxn_ids_serine.txt', 
                     header = F, sep = '\t', quote = '')

flux = cbind(flux_denovo, flux_plasma, flux_cortex)
no_rxns = nrow(rxn_ids)
flux = t(flux)
head(rownames(flux))
dim(flux)
colnames(flux) = paste(rep('v', no_rxns), seq(1, no_rxns), sep = '')

flux_names = paste(rep('v', no_rxns), seq(1, no_rxns), sep = '')
flux = as.data.frame(flux)
flux[['sim']] = factor(rownames(flux))

flux_long = gather(flux, fluxNo, value, v1:v21, factor_key = T)

net_rxn = c('v1', 'v3')
net_rxn_ratio = c('v2', 'v4')
irrev_rxn = paste(rep('v', no_rxns - 7 + 1), seq(7, no_rxns), sep = '')
exchange_rxn = paste(rep('v', no_rxns - 16 + 1), seq(16, no_rxns), sep = '')

irrev_rxn = irrev_rxn[!irrev_rxn %in% exchange_rxn]

hex_codes = c()
for (x in flux_names){
  
  if (x %in% net_rxn){
    hex_codes = c(hex_codes, '#C70039')
  } else if (x %in% net_rxn_ratio){
    hex_codes = c(hex_codes, '#32CD32')
  } else if (x %in% irrev_rxn){
    hex_codes = c(hex_codes, '#FFC300')
  } else {
    hex_codes = c(hex_codes, '#355C7D')
  }
}

flux_long_selected_1 = flux_long[flux_long$fluxNo %in% c(net_rxn, irrev_rxn, 
                                                         exchange_rxn),]
flux_long_selected_2 = flux_long[flux_long$fluxNo %in% net_rxn_ratio,]

for (i in unique(flux_long_selected_1$fluxNo)){
  flux_long_selected_1[flux_long_selected_1$fluxNo == i, 
                       'reactions'] = bounds[bounds$flux_no == i, 'reactions']
}

# Visualize flux distribution for net reversible, irreversible, and sink fluxes
p = flux_long_selected_1 %>%
  ggplot() + 
  geom_violinhalf(aes(x = fluxNo, y = value, fill = fluxNo), 
                  trim = F, lwd = 0.3, scale = 'width') +
  
  geom_half_boxplot(aes(x = fluxNo, y = value),
                    width = 0.45, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  
  theme_classic()+
  labs(
    x    = "Flux Index",
    y    = "Simulated Flux Distributions",
    fill = "Flux Index"
  ) +
  scale_fill_manual(values = unname(hex_codes[unique(flux_long_selected_1$fluxNo)])) +
  coord_flip() + 
  theme(legend.position = 'none')

ggsave(plot = p, filename = paste('flux_distribution', 
                                  'glioma_50_plasma_denovo_cortex.pdf', 
                                  sep = '_'),
       device = "pdf", height = 4, 
       width = 3, units = "in")

# Visualize distribution of exchange factor of reversible reactions
p = flux_long_selected_2 %>%
  ggplot() + 
  geom_violinhalf(aes(x = fluxNo, y = value, fill = fluxNo), 
                  trim = T, lwd = 0.3, scale = 'width') +
  
  geom_half_boxplot(aes(x = fluxNo, y = value),
                    width = 0.45, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 

  theme_classic()+
  labs(
    x    = "Flux Index",
    y    = "Simulated Flux Distributions",
    fill = "Flux Index"
  ) +
  scale_fill_manual(values = unname(hex_codes[unique(flux_long_selected_1$fluxNo)])) +
  coord_flip() + 
  theme(legend.position = 'none')

ggsave(plot = p, filename = paste('flux_distribution', 
                                  'glioma_50_plasma_denovo_cortex_2.pdf', 
                                  sep = '_'),
       device = "pdf", height = 1, 
       width = 1.5, units = "in")


# Purine flux distribution ----------------------------------------------------

# Read simulated fluxes
flux_denovo = read.delim(file = './sim_flux/flux_glioma_denovo_gmp.txt', 
                         header = F, sep = '\t', quote = '')

rxn_ids = read.delim(file = './sim_flux/rxn_ids_purine.txt', 
                     header = F, sep = '\t', quote = '')

no_rxns = nrow(rxn_ids)
flux_denovo = t(flux_denovo)
head(rownames(flux_denovo))
dim(flux_denovo)
colnames(flux_denovo) = paste(rep('v', no_rxns), seq(1, no_rxns), sep = '')

flux_names = paste(rep('v', no_rxns), seq(1, no_rxns), sep = '')
flux_denovo = as.data.frame(flux_denovo)
flux_denovo[['sim']] = factor(rownames(flux_denovo))

flux_long = tidyr::gather(flux_denovo, fluxNo, value, v1:v22, factor_key = T)

net_rxn = c('v1')
net_rxn_ratio = c('v2')
irrev_rxn = paste(rep('v', no_rxns - 3 + 1), seq(3, no_rxns), sep = '')
exchange_rxn = c('v4', 'v5', 'v11', 'v14', 'v19', 'v22')

irrev_rxn = irrev_rxn[!irrev_rxn %in% exchange_rxn]

hex_codes = c()
for (x in flux_names){
  
  if (x %in% net_rxn){
    hex_codes = c(hex_codes, '#C70039')
  } else if (x %in% net_rxn_ratio){
    hex_codes = c(hex_codes, '#32CD32')
  } else if (x %in% irrev_rxn){
    hex_codes = c(hex_codes, '#FFC300')
  } else {
    hex_codes = c(hex_codes, '#355C7D')
  }
}

flux_long_selected_1 = flux_long[flux_long$fluxNo %in% c(net_rxn, irrev_rxn, 
                                                         exchange_rxn),]
flux_long_selected_2 = flux_long[flux_long$fluxNo %in% net_rxn_ratio,]

# Visualize flux distribution for net reversible, irreversible, and sink fluxes
p = flux_long_selected_1 %>%
  ggplot() + 
  geom_violinhalf(aes(x = fluxNo, y = value, fill = fluxNo), 
                  trim = F, lwd = 0.3, scale = 'width') +
  
  geom_half_boxplot(aes(x = fluxNo, y = value),
                    width = 0.45, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  
  theme_classic()+
  labs(
    x    = "Flux Index",
    y    = "Random Flux Distributions",
    fill = "Flux Index"
  ) +
  scale_fill_manual(values = unname(hex_codes[unique(flux_long_selected_1$fluxNo)])) +
  coord_flip() + 
  theme(legend.position = 'none')

ggsave(plot = p, filename = paste('flux_distribution',
                                  'gbm_purine_denovo.pdf', sep = '_'),
       device = "pdf", height = 4, 
       width = 3, units = "in")

# Visualize distribution of exchange factor of reversible reactions
p = flux_long_selected_2 %>%
  ggplot() + 
  geom_violinhalf(aes(x = fluxNo, y = value, fill = fluxNo), 
                  trim = T, lwd = 0.3, scale = 'width') +
  
  geom_half_boxplot(aes(x = fluxNo, y = value),
                    width = 0.45, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  theme_classic()+
  labs(
    x    = "Flux Index",
    y    = "Simulated Flux Distributions",
    fill = "Flux Index"
  ) +
  scale_fill_manual(values = unname(hex_codes[unique(flux_long_selected_1$fluxNo)])) +
  coord_flip() + 
  theme(legend.position = 'none')

ggsave(plot = p, filename = paste('flux_distribution', 'gbm_purine_denovo_2.pdf', 
                                  sep = '_'),
       device = "pdf", height = 1,
       width = 1.5, units = "in")



# Plots for purine simulated MIDs --------------------------------------------

## Read purine simulated MIDs ----
rm(list = ls())

mid = read.csv(file = paste('../metabolic_CNN/data/sim_data/', 
                            'simulated_data_gmp_denovo_glioma.csv', 
                            sep = ''), header = T, quote = '')

patient_mids = read.delim(file = 'patient_MIDs_cortex_glioma_purine.txt', 
                          header = T, sep = '\t', quote = '')

patient_mids$INFUSION_TIME = patient_mids$INFUSION_TIME/60 # hour

max(patient_mids$INFUSION_TIME) # 3.66 hr
min(patient_mids$INFUSION_TIME) # 2.2 hr

patient_time_points = unique(round(patient_mids$INFUSION_TIME, digits = 1))
patient_time_points

patient_mid_rm = subset(patient_mids, 
                        select = -c(R5P0, R5P1, R5P2, R5P3, R5P4, R5P5,
                                    IMP0, GMP0, INO0, GUO0, GDP0, AMP0))

head(patient_mid_rm)

patient_mid_rm$PATIENT = str_sub(patient_mid_rm$PATIENT, start = -1)
patient_mid_rm$SITENO = str_sub(patient_mid_rm$SITE, start = -1)
patient_mid_rm$SITE = str_sub(patient_mid_rm$SITE, start = 1, end = -2)

patient_mid_rm$INFUSION_TIME = round(patient_mid_rm$INFUSION_TIME, digits = 1)

patient_mid_long = pivot_longer(patient_mid_rm, cols = 5:ncol(patient_mid_rm)-1)

head(patient_mid_long)

patient_mid_long[['PATIENTSITE']] = paste(patient_mid_long$PATIENT, 
                                          patient_mid_long$SITE, sep = '')

sim_mid = mid[mid$time %in% patient_time_points, ]
head(sim_mid)

sim_mid[, seq(3, ncol(sim_mid))] = sim_mid[, seq(3, ncol(sim_mid))] * 100

sim_mid_rmM0 = subset(sim_mid, select = -c(GLY0, GLY1, GLY2,
                                           MTHF0, MTHF1,
                                           IMP0, GMP0, INO0, GUO0, GDP0, AMP0))
head(sim_mid_rmM0)
ncol(sim_mid_rmM0)

sim_mid_long = pivot_longer(sim_mid_rmM0, cols = 3:ncol(sim_mid_rmM0))
sim_mid_long[['Mi']] = str_sub(sim_mid_long$name, start = -1)
colors = c("#FFC300", "#DA8A67", "#E41B17", "#7D0552", "#808080")

for (i in seq(1, 5)){
  sim_mid_long[sim_mid_long$Mi == i, 'colors'] = colors[i]  
}


gc()
## Visualize distribution of MIDs at patient time points ----
for (i in patient_time_points){
  sim_mid_t = sim_mid_long[sim_mid_long$time == i, ]
  patient_points = patient_mid_long[patient_mid_long$INFUSION_TIME == i, ]
  p = ggplot() + 
    
    geom_violinhalf(data = sim_mid_t, aes(x = name, y = value, fill = Mi), 
                    trim = T, lwd = 0.3, scale = 'width') +
    
    geom_half_boxplot(data = sim_mid_t, aes(x = name, y = value),
                      width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
    geom_point(data = patient_points, aes(x = name, y = value, colour = PATIENTSITE),
               shape = 1, size = 0.1) + 
    
    theme_classic()+
    labs(
      x    = "Metabolite MID",
      y    = "Simulated MIDs",
      fill = "M+i"
    ) +
    scale_fill_manual(values = colors) +
    coord_flip() + 
    theme(legend.key.size = unit(0.05, 'in'), legend.title = element_text(size=7))
  
  ggsave(plot = p, filename = paste('sim_MID', i, '.pdf', sep = '_'),
         device = "pdf", height = 4, 
         width = 4, units = "in")
  
}

## Visualize distribution of MIDs all timepoints ----
p = ggplot() + 
  
  geom_violinhalf(data = sim_mid_long, aes(x = name, y = value, fill = Mi), 
                  trim = T, lwd = 0.3, scale = 'width') +
  
  geom_half_boxplot(data = sim_mid_long, aes(x = name, y = value),
                    width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_point(data = patient_mid_long, aes(x = name, y = value, colour = PATIENTSITE),
             shape = 1, size = 0.1) + 
  theme_classic()+
  labs(
    x    = "Metabolite MID",
    y    = "Simulated MIDs",
    fill = "M+i"
  ) +
  scale_fill_manual(values = colors) +
  coord_flip() + 
  theme(legend.key.size = unit(0.05, 'in'), legend.title = element_text(size=7))

ggsave(plot = p, 
       filename = 'sim_MID_all_patients_trim_purine_denovo_gbm.pdf', 
       device = "pdf", height = 5.5,
       width = 4, units = "in")



# tsne plots for purine MID simulation ----------------------------------------

## Read simulation parameters including serine and R5P ----

param_mid = read.csv(file =  paste('../metabolic_CNN/data/sim_data/', 
                             'simulation_parameters_gmp_denovo_glioma.csv', 
                             sep = ''), header = T, quote = '')

param_mid[['index']] = seq(1, dim(param_mid)[1])

add_mid <- function(mid_df, added_mid, param_df){
  mid_df[[added_mid]] = param_df[match(mid_df$index, param_df$index), added_mid]
  return(mid_df)
}

## Combine plasma serine and R5P MIDs from parameters to mid dataframe ----
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

mid = mid[, !(colnames(mid) %in% c('MTHF0', 'MTHF1', 'GLY0', 'GLY1', 'GLY2'))]

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
  
  mid_v = rbind(mid_mc_ps, mid_v)
  
}

mid_v = data.matrix(mid_v)
mid_v = t(mid_v)

com_mid = CreateSeuratObject(counts = mid_v, 
                             project = "combined_sim_patient_mid_purine")

com_mid = FindVariableFeatures(com_mid, selection.method = "vst")

com_mid = ScaleData(com_mid, features = rownames(com_mid))

com_mid = RunPCA(com_mid, 
                 features = VariableFeatures(object = com_mid), npcs = 10)

ElbowPlot(com_mid)

com_mid = FindNeighbors(com_mid, dims = 1:5)

com_mid = FindClusters(com_mid, resolution = 0.1)

com_mid = RunTSNE(com_mid, dims = 1:5)

my_colors = c('P8N' = "#5A5156",
              'P7N' = "#C4451C",
              "P7E" = "#822E1C",
              "P6N" = "#1C7F93",
              "P6E" = "#3B00FB",
              "P5N" = "#D85FF7",
              "P5E" = "#FEAF16",
              "P4N" = "#90AD1C",
              "P4E" = "#1C8356",
              "P3N" = "#AA0DFE",
              "P3E" = "#DEA0FD",
              "P2N" = "#2ED9FF",
              "P2E" = "#3283FE",
              "P1N" = "#B00068",
              "P1E" = "#F6222E",
              'sim' = "#E4E1E3")

PS = c(rep('P8N', 100), rep('P7N', 100), rep('P7E', 100), 
       rep('P6N', 100), rep('P6E', 100), rep('P5N', 100), rep('P5E', 100),
       rep('P4N', 100), rep('P4E', 100), rep('P3N', 100), rep('P3E', 100),
       rep('P2N', 100), rep('P2E', 100), rep('P1N', 100), rep('P1E', 100),
       rep('sim', 1050000))

com_mid@meta.data[['ps_sim']] = PS

saveRDS(com_mid, 'sim_ps_glioma_50_purine.rds')

pdf(file = 'dimplot_tsne_sim_purine.pdf', height = 10, width = 10)
DimPlot(com_mid, reduction = "tsne", group.by = 'ps_sim', 
        cells = colnames(com_mid)[1501:1051500], cols = my_colors)
dev.off()

pdf(file = 'dimplot_tsne_patient_site_purine.pdf', height = 10, width = 10)
DimPlot(com_mid, reduction = "tsne", group.by = 'ps_sim', 
        cells = colnames(com_mid)[1:1500], cols = my_colors)
dev.off()

