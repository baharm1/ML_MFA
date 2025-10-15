library(tidyr)
library(RColorBrewer)
library(see)
library(ggplot2)
library(dplyr)
library(gghalves)
library(stringr)

rm(list = ls())
# simulated_data_v9.csv with changed flux bounds paramater_bounds_v3
# simulated_data_v8.csv with Anjali's flux bounds parameter_bounds_v1
version = 'v9_serine'
mid = read.csv(file = paste('../matlab_files_v3_serine/serine_v9_cortex_50/simulated_data_', version, '.csv', 
                            sep = ''), header = T, quote = '')
patient_mids = read.delim(file = '../patient_input_data/patient_MIDs_serine_added.txt', 
                          header = T, sep = '\t', quote = '')

patient_mids$INFUSION_TIME = patient_mids$INFUSION_TIME/60 # hour

max(patient_mids$INFUSION_TIME) # 3.66 hr
min(patient_mids$INFUSION_TIME) # 2.2 hr

patient_time_points = unique(round(patient_mids$INFUSION_TIME, digits = 1))
patient_time_points

#patient_mids[['psite']] = paste(patient_mids$PATIENT, patient_mids$SITE,sep = '')

patient_mid_rm = subset(patient_mids, select = -c(IMP0, IMP1, IMP2, IMP3, IMP4, IMP5,
                                                  GDP0, GDP1, GDP2, GDP3, GDP4, GDP5,
                                                  GUO0, GUO1, GUO2, GUO3, GUO4, GUO5,
                                                  INO0, INO1, INO2, INO3, INO4, INO5,
                                                  AMP0, AMP1, AMP2, AMP3, AMP4, AMP5,
                                                  GMP0, GMP1, GMP2, GMP3, GMP4, GMP5,
                                                  PG0, PGa0, SERa0, SER0, 
                                                  GLY0, R5P0))


head(patient_mid_rm)

patient_mid_rm$PATIENT = str_sub(patient_mid_rm$PATIENT, start = -1)
patient_mid_rm$SITENO = str_sub(patient_mid_rm$SITE, start = -1)
patient_mid_rm$SITE = ifelse(str_sub(patient_mid_rm$SITE, start = 1, end = -2) == 'E', 1, 0)

patient_mid_rm$INFUSION_TIME = round(patient_mid_rm$INFUSION_TIME, digits = 1)

patient_mid_long = pivot_longer(patient_mid_rm, cols = 5:ncol(patient_mid_rm)-1)
#patient_mid_rm = patient_mid_rm[, -c(ncol(patient_mid_rm)-1)]
head(patient_mid_long)

patient_mid_long[['PATIENTSITE']] = paste(patient_mid_long$PATIENT, patient_mid_long$SITE, sep = '')
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

sim_mid_rmM0 = subset(sim_mid, select = -c(GLY0, SER0, SERa0, CTHF0, CTHF1, PG0,
                                           R5P0, AMP0, GDP0, GMP0, GUO0, 
                                           IMP0, INO0))
head(sim_mid_rmM0)
ncol(sim_mid_rmM0)

sim_mid_long = pivot_longer(sim_mid_rmM0, cols = 3:ncol(sim_mid_rmM0))
sim_mid_long[['Mi']] = str_sub(sim_mid_long$name, start = -1)
colors = c("#FFC300", "#DA8A67", "#E41B17", "#7D0552", "#2F0909")

for (i in seq(1, 5)){
  sim_mid_long[sim_mid_long$Mi == i, 'colors'] = colors[i]  
}


gc()

for (i in patient_time_points){
  sim_mid_t = sim_mid_long[sim_mid_long$time == i, ]
  patient_points = patient_mid_long[patient_mid_long$INFUSION_TIME == i, ]
  p = ggplot() + 
    #ggdist::stat_halfeye(adjust = 0.5, justification = -0.3, .width = 0, point_colour = NA, scale = 0.5) +
    geom_violinhalf(data = sim_mid_t, aes(x = name, y = value, fill = Mi), 
                    trim = F, lwd = 0.3, scale = 'width') +
    
    #geom_boxplot(aes(x = fluxNo, y = value),
    #            width = 0.2, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
    geom_half_boxplot(data = sim_mid_t, aes(x = name, y = value),
                      width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
    geom_point(data = patient_points, aes(x = name, y = value, colour = PATIENTSITE),
               shape = 1, size = 0.1) + 
    #geom_point(aes(x = 'SER1', y = 0.1), shape = 2, size = 0.1) + 
    
    theme_classic()+
    labs(
      #title = "Raincloud Plot",
      #subtitle = "Showing the Multi-Modal Distributions of Temperatures by Month",
      x    = "Metabolite MID",
      y    = "Simulated MIDs",
      fill = "M+i"
    ) +
    #ggbreak::scale_y_break(c(20, 40)) + # 33 v9, 40 v8
    scale_fill_manual(values = colors) +
    #ylim(-0.5, 20) +
    coord_flip() + 
    theme(legend.key.size = unit(0.05, 'in'), legend.title = element_text(size=7))
  
    #theme(legend.position = 'none')
  
  ggsave(plot = p, filename = paste(version, 'sim_MID', i, '.pdf', sep = '_'),
         device = "pdf", height = 7, #5.5
         width = 5, units = "in")#3.2
  
}

# plots for all time points and all patients ----
rm(list = ls())
# simulated_data_v9.csv with changed flux bounds paramater_bounds_v3
# simulated_data_v8.csv with Anjali's flux bounds parameter_bounds_v1
version = 'v29'
mid = read.csv(file = paste('../matlab_files_v2/simulated_data_', version, '.csv', 
                            sep = ''), header = T, quote = '')
patient_mids = read.delim(file = '../patient_input_data/patient_MIDs_serine_added.txt', 
                          header = T, sep = '\t', quote = '')

timepoints_ML = seq(2, 4, 0.1)
sim_mid = mid[mid$time %in% timepoints_ML, ]

patient_mids$INFUSION_TIME = patient_mids$INFUSION_TIME/60 # hour

max(patient_mids$INFUSION_TIME) # 3.66 hr
min(patient_mids$INFUSION_TIME) # 2.2 hr

patient_time_points = unique(round(patient_mids$INFUSION_TIME, digits = 1))
patient_time_points

#patient_mids[['psite']] = paste(patient_mids$PATIENT, patient_mids$SITE,sep = '')

patient_mid_rm = subset(patient_mids, select = -c(IMP0, GDP0, GUO0, INO0, AMP0, 
                                                  GMP0, PG0, PGa0, PGa1, PGa2, PGa3,
                                                  SERa0, SER0, GLY0, R5P0))

head(patient_mid_rm)

patient_mid_rm$PATIENT = str_sub(patient_mid_rm$PATIENT, start = -1)
patient_mid_rm$SITENO = str_sub(patient_mid_rm$SITE, start = -1)
patient_mid_rm$SITE = ifelse(str_sub(patient_mid_rm$SITE, start = 1, end = -2) == 'E', 1, 0)

patient_mid_rm$INFUSION_TIME = round(patient_mid_rm$INFUSION_TIME, digits = 1)

patient_mid_long = pivot_longer(patient_mid_rm, cols = 5:ncol(patient_mid_rm)-1)
#patient_mid_rm = patient_mid_rm[, -c(ncol(patient_mid_rm)-1)]
head(patient_mid_long)

patient_mid_long[['PATIENTSITE']] = paste(patient_mid_long$PATIENT, patient_mid_long$SITE, sep = '')
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
  

#sim_mid = mid[mid$time %in% patient_time_points, ]
#head(sim_mid)

sim_mid[, seq(3, ncol(sim_mid))] = sim_mid[, seq(3, ncol(sim_mid))] * 100

sim_mid_rmM0 = subset(sim_mid, select = -c(GLY0, SER0, SERa0, CTHF0, CTHF1, PG0,
                                           R5P0, AMP0, GDP0, GMP0, GUO0, 
                                           IMP0, INO0))
head(sim_mid_rmM0)
ncol(sim_mid_rmM0)

sim_mid_long = pivot_longer(sim_mid_rmM0, cols = 3:ncol(sim_mid_rmM0))
sim_mid_long[['Mi']] = str_sub(sim_mid_long$name, start = -1)
colors = c("#FFC300", "#DA8A67", "#E41B17", "#7D0552", "#2F0909")

for (i in seq(1, 5)){
  sim_mid_long[sim_mid_long$Mi == i, 'colors'] = colors[i]  
}


gc()

p = ggplot() + 
  #ggdist::stat_halfeye(adjust = 0.5, justification = -0.3, .width = 0, point_colour = NA, scale = 0.5) +
  geom_violinhalf(data = sim_mid_long, aes(x = name, y = value, fill = Mi), 
                  trim = T, lwd = 0.3, scale = 'width') +
  
  #geom_boxplot(aes(x = fluxNo, y = value),
  #            width = 0.2, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_half_boxplot(data = sim_mid_long, aes(x = name, y = value),
                    width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_point(data = patient_mid_long, aes(x = name, y = value, colour = PATIENTSITE),
             shape = 1, size = 0.1) + 
  #geom_point(aes(x = 'SER1', y = 0.1), shape = 2, size = 0.1) + 
  
  theme_classic()+
  labs(
    #title = "Raincloud Plot",
    #subtitle = "Showing the Multi-Modal Distributions of Temperatures by Month",
    x    = "Metabolite MID",
    y    = "Simulated MIDs",
    fill = "M+i"
  ) +
  #ggbreak::scale_y_break(c(20, 40)) + # 33 v9, 40 v8
  scale_fill_manual(values = colors) +
  #ylim(-0.5, 20) +
  coord_flip() + 
  theme(legend.key.size = unit(0.05, 'in'), legend.title = element_text(size=7))

#theme(legend.position = 'none')

ggsave(plot = p, filename = paste(version, 'sim_MID_all_patients_trim.pdf', sep = '_'),
       device = "pdf", height = 7, #5.5
       width = 5, units = "in")#3.2



# plots for serine simulation ----
library(tidyr)
library(RColorBrewer)
library(see)
library(ggplot2)
library(dplyr)
library(gghalves)
library(stringr)

rm(list = ls())
# simulated_data_v9.csv with changed flux bounds paramater_bounds_v3
# simulated_data_v8.csv with Anjali's flux bounds parameter_bounds_v1
version = 'v9_serine'
mid_plasma = read.csv(file = paste('../matlab_files_v3_serine/serine_v9_glioma_50_plasma/simulated_data_', version, '.csv', 
                            sep = ''), header = T, quote = '')
mid_denovo = read.csv(file = paste('../matlab_files_v3_serine/serine_v9_glioma_50_denovo/simulated_data_', version, '.csv', 
                                   sep = ''), header = T, quote = '')
mid_cortex = read.csv(file = paste('../matlab_files_v3_serine/serine_v9_cortex_50/simulated_data_', version, '.csv', 
                                   sep = ''), header = T, quote = '')

mid = rbind(mid_denovo, mid_plasma, mid_cortex)

patient_mids = read.delim(file = '../patient_input_data/patient_MIDs_serine_added_cortex_gbm.txt', 
                          header = T, sep = '\t', quote = '')

patient_mids$INFUSION_TIME = patient_mids$INFUSION_TIME/60 # hour

max(patient_mids$INFUSION_TIME) # 3.66 hr
min(patient_mids$INFUSION_TIME) # 2.2 hr

patient_time_points = unique(round(patient_mids$INFUSION_TIME, digits = 1))
patient_time_points

#patient_mids[['psite']] = paste(patient_mids$PATIENT, patient_mids$SITE,sep = '')

patient_mid_rm = subset(patient_mids, select = -c(PG0, PGc0, SER0, SERc0, GLY0, GLYc0))

head(patient_mid_rm)

patient_mid_rm$PATIENT = str_sub(patient_mid_rm$PATIENT, start = -1)
patient_mid_rm$SITENO = str_sub(patient_mid_rm$SITE, start = -1)
patient_mid_rm$SITE = ifelse(str_sub(patient_mid_rm$SITE, start = 1, end = -2) == 'E', 1, 0)

patient_mid_rm$INFUSION_TIME = round(patient_mid_rm$INFUSION_TIME, digits = 1)

patient_mid_long = pivot_longer(patient_mid_rm, cols = 5:ncol(patient_mid_rm)-1)
#patient_mid_rm = patient_mid_rm[, -c(ncol(patient_mid_rm)-1)]
head(patient_mid_long)

patient_mid_long[['PATIENTSITE']] = paste(patient_mid_long$PATIENT, patient_mid_long$SITE, sep = '')
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
                                           MTHF0, MTHF1, MTHFc0, MTHFc1, PG0, PGc0))
head(sim_mid_rmM0)
ncol(sim_mid_rmM0)

sim_mid_long = pivot_longer(sim_mid_rmM0, cols = 3:ncol(sim_mid_rmM0))
sim_mid_long[['Mi']] = str_sub(sim_mid_long$name, start = -1)
colors = c("#FFC300", "#DA8A67", "#E41B17", "#7D0552", "#2F0909")

for (i in seq(1, 5)){
sim_mid_long[sim_mid_long$Mi == i, 'colors'] = colors[i]  
}


gc()

for (i in patient_time_points){
sim_mid_t = sim_mid_long[sim_mid_long$time == i, ]
patient_points = patient_mid_long[patient_mid_long$INFUSION_TIME == i, ]
p = ggplot() + 
  #ggdist::stat_halfeye(adjust = 0.5, justification = -0.3, .width = 0, point_colour = NA, scale = 0.5) +
  geom_violinhalf(data = sim_mid_t, aes(x = name, y = value, fill = Mi), 
                  trim = T, lwd = 0.3, scale = 'width') +
  
  #geom_boxplot(aes(x = fluxNo, y = value),
  #            width = 0.2, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_half_boxplot(data = sim_mid_t, aes(x = name, y = value),
                    width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_point(data = patient_points, aes(x = name, y = value, colour = PATIENTSITE),
             shape = 1, size = 0.1) + 
  #geom_point(aes(x = 'SER1', y = 0.1), shape = 2, size = 0.1) + 
  
  theme_classic()+
  labs(
    #title = "Raincloud Plot",
    #subtitle = "Showing the Multi-Modal Distributions of Temperatures by Month",
    x    = "Metabolite MID",
    y    = "Simulated MIDs",
    fill = "M+i"
  ) +
  #ggbreak::scale_y_break(c(20, 40)) + # 33 v9, 40 v8
  scale_fill_manual(values = colors) +
  #ylim(-0.5, 30) +
  coord_flip() + 
  theme(legend.key.size = unit(0.05, 'in'), legend.title = element_text(size=7))

#theme(legend.position = 'none')

ggsave(plot = p, filename = paste(version, 'sim_MID', i, '1.pdf', sep = '_'),
       device = "pdf", height = 4, #5.5
       width = 4, units = "in")#3.2

}


p = ggplot() + 
  #ggdist::stat_halfeye(adjust = 0.5, justification = -0.3, .width = 0, point_colour = NA, scale = 0.5) +
  geom_violinhalf(data = sim_mid_long, aes(x = name, y = value, fill = Mi), 
                  trim = T, lwd = 0.3, scale = 'width') +
  
  #geom_boxplot(aes(x = fluxNo, y = value),
  #            width = 0.2, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_half_boxplot(data = sim_mid_long, aes(x = name, y = value),
                    width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_point(data = patient_mid_long, aes(x = name, y = value, colour = PATIENTSITE),
             shape = 1, size = 0.1) + 
  #geom_point(aes(x = 'SER1', y = 0.1), shape = 2, size = 0.1) + 
  
  theme_classic()+
  labs(
    #title = "Raincloud Plot",
    #subtitle = "Showing the Multi-Modal Distributions of Temperatures by Month",
    x    = "Metabolite MID",
    y    = "Simulated MIDs",
    fill = "M+i"
  ) +
  #ggbreak::scale_y_break(c(20, 40)) + # 33 v9, 40 v8
  scale_fill_manual(values = colors) +
  #ylim(-0.5, 20) +
  coord_flip() + 
  theme(legend.key.size = unit(0.05, 'in'), legend.title = element_text(size=7))

#theme(legend.position = 'none')

ggsave(plot = p, filename = paste(version, 'sim_MID_all_patients_trim_all_glioma_denovo_plasma_cortex_50.pdf', sep = '_'),
       device = "pdf", height = 4, #5.5
       width = 4, units = "in")#3.2

# tsne plots for serine MID simulation ----
# simulation MIDs
version = 'v9_serine'
mid_plasma = read.csv(file = paste('../matlab_files_v3_serine/serine_v9_glioma_50_plasma/simulated_data_', version, '.csv', 
                                   sep = ''), header = T, quote = '')
mid_denovo = read.csv(file = paste('../matlab_files_v3_serine/serine_v9_glioma_50_denovo/simulated_data_', version, '.csv', 
                                   sep = ''), header = T, quote = '')
mid_cortex = read.csv(file = paste('../matlab_files_v3_serine/serine_v9_cortex_50/simulated_data_', version, '.csv', 
                                   sep = ''), header = T, quote = '')

# simulation parameters including plasma serine
param_plasma = read.csv(file = paste('../matlab_files_v3_serine/serine_v9_glioma_50_plasma/simulation_parameters_', version, '.csv', 
                                     sep = ''), header = T, quote = '')
param_denovo = read.csv(file = paste('../matlab_files_v3_serine/serine_v9_glioma_50_denovo/simulation_parameters_', version, '.csv', 
                                   sep = ''), header = T, quote = '')
param_cortex = read.csv(file = paste('../matlab_files_v3_serine/serine_v9_cortex_50/simulation_parameters_', version, '.csv', 
                                   sep = ''), header = T, quote = '')

param_plasma[['index']] = seq(1, dim(param_plasma)[1])
param_denovo[['index']] = seq(1, dim(param_denovo)[1])
param_cortex[['index']] = seq(1, dim(param_cortex)[1])

add_mid <- function(mid_df, added_mid, param_df){
  mid_df[[added_mid]] = param_df[match(mid_df$index, param_df$index), added_mid]
  return(mid_df)
}

# combine plasma serine MIDs from parameters to mid dataframe
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

# patient MIDs
patient_dir = '../patient_MFA/patient_serine_MFA/Output files/serine_MFA_second_try_scflux'
patient_files = list.dirs(path = patient_dir, recursive = F)
mid_names = read.delim('../patient_MFA/patient_serine_MFA/mid_name.txt', sep = '\t', header = F)

# combine simulated mids with patient mids
for (ps in 1:length(patient_files)){
  
  mid_mc_ps = read.delim(paste(patient_files[ps], '/mid_mc.txt', sep = ''), 
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

pdf(file = 'dimplot_tsne_plasma_sim_pat_com_PS_3.pdf', height = 10, width = 10)
DimPlot_scCustom(com_plasma, reduction = "tsne", group.by = 'ps_sim', 
        cells = colnames(com_plasma)[1501:1051500], colors_use = my_colors)
dev.off()

pdf(file = 'dimplot_tsne_plasma_sim_pat_com_PS_6.pdf', height = 10, width = 10)
DimPlot(com_plasma, reduction = "tsne", group.by = 'ps_sim', 
                 cells = colnames(com_plasma)[1501:1051500], cols = my_colors)
dev.off()

PS = c(rep('P8N', 100), rep('P7N', 100), rep('P7E', 100), 
       rep('P6N', 100), rep('P6E', 100), rep('P5N', 100), rep('P5E', 100),
       rep('P4N', 100), rep('P4E', 100), rep('P3N', 100), rep('P3E', 100),
       rep('P2N', 100), rep('P2E', 100), rep('P1N', 100), rep('P1E', 100),
       rep('sim', 1050000))

com_plasma@meta.data[['ps_sim']] = PS

saveRDS(com_plasma, 'sim_glioma_50_plasma_ps.rds')

# serine flux distribution ----
flux_plasma = read.delim(file = paste('../matlab_files_v3_serine/serine_v9_glioma_50_plasma/flux_', 
                               version, '.txt', sep = ''),
                  header = F, sep = '\t', quote = '')

rxn_ids = read.delim(file = paste('../matlab_files_v3_serine/serine_v9_glioma_50_plasma/rxn_ids_', 
                                  version, '.txt', sep = ''),
                     header = F, sep = '\t', quote = '')

flux_denovo = read.delim(file = paste('../matlab_files_v3_serine/serine_v9_glioma_50_denovo/flux_', 
                                      version, '.txt', sep = ''),
                         header = F, sep = '\t', quote = '')

flux_cortex = read.delim(file = paste('../matlab_files_v3_serine/serine_v9_cortex_50/flux_', 
                                      version, '.txt', sep = ''),
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

#exchange_rxn = c('v32', 'v33', 'v34', 'v35', 'v36', 'v37', 'v38', 
#'v44', 'v47', 'v52', 'v53')
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
  flux_long_selected_1[flux_long_selected_1$fluxNo == i, 'reactions'] = bounds[bounds$flux_no == i, 'reactions']
}


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
  #ggbreak::scale_y_break(c(1150, 8400)) +
  scale_fill_manual(values = unname(hex_codes[unique(flux_long_selected_1$fluxNo)])) +
  coord_flip() + 
  theme(legend.position = 'none')

ggsave(plot = p, filename = paste('flux_distribution', version, 'glioma_50_plasma_denovo_cortex.pdf', sep = '_'),
       device = "pdf", height = 4, #5.5
       width = 3, units = "in")#3.2

ggsave(plot = p, filename = paste('flux_distribution', version, 'glioma_50_plasma_denovo_cortex.jpg', sep = '_'),
       device = "jpeg", height = 6, #5.5
       width = 3.2, units = "in", dpi = 1000)#3.2

p = flux_long_selected_2 %>%
  ggplot() + 
  geom_violinhalf(aes(x = fluxNo, y = value, fill = fluxNo), 
                  trim = T, lwd = 0.3, scale = 'width') +
  
  geom_half_boxplot(aes(x = fluxNo, y = value),
                    width = 0.45, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  #geom_point(aes(x = fluxNo, y = lb), shape = 6, size = 0.1) + 
  #geom_point(aes(x = fluxNo, y = ub), shape = 2, size = 0.1) + 
  
  theme_classic()+
  labs(
    x    = "Flux Index",
    y    = "Random Flux Distributions",
    fill = "Flux Index"
  ) +
  #ggbreak::scale_y_break(c(1150, 8400)) +
  scale_fill_manual(values = unname(hex_codes[unique(flux_long_selected_1$fluxNo)])) +
  coord_flip() + 
  theme(legend.position = 'none')

ggsave(plot = p, filename = paste('flux_distribution', version, 'glioma_50_plasma_denovo_cortex_2.pdf', sep = '_'),
       device = "pdf", height = 1, #5.5
       width = 1.5, units = "in")#3.2


# purine flux distribution ----
version = 'v5'
flux_denovo = read.delim(file = paste('../matlab_files_purine/flux_', 
                                      version, '.txt', sep = ''),
                         header = F, sep = '\t', quote = '')

rxn_ids = read.delim(file = paste('../matlab_files_purine/rxn_ids_', 
                                  version, '.txt', sep = ''),
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

#exchange_rxn = c('v32', 'v33', 'v34', 'v35', 'v36', 'v37', 'v38', 
#'v44', 'v47', 'v52', 'v53')
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
  #ggbreak::scale_y_break(c(1150, 8400)) +
  scale_fill_manual(values = unname(hex_codes[unique(flux_long_selected_1$fluxNo)])) +
  coord_flip() + 
  theme(legend.position = 'none')

ggsave(plot = p, filename = paste('flux_distribution', version, 'gbm_purine_denovo.pdf', sep = '_'),
       device = "pdf", height = 4, #5.5
       width = 3, units = "in")#3.2


p = flux_long_selected_2 %>%
  ggplot() + 
  geom_violinhalf(aes(x = fluxNo, y = value, fill = fluxNo), 
                  trim = T, lwd = 0.3, scale = 'width') +
  
  geom_half_boxplot(aes(x = fluxNo, y = value),
                    width = 0.45, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  #geom_point(aes(x = fluxNo, y = lb), shape = 6, size = 0.1) + 
  #geom_point(aes(x = fluxNo, y = ub), shape = 2, size = 0.1) + 
  
  theme_classic()+
  labs(
    x    = "Flux Index",
    y    = "Random Flux Distributions",
    fill = "Flux Index"
  ) +
  #ggbreak::scale_y_break(c(1150, 8400)) +
  scale_fill_manual(values = unname(hex_codes[unique(flux_long_selected_1$fluxNo)])) +
  coord_flip() + 
  theme(legend.position = 'none')

ggsave(plot = p, filename = paste('flux_distribution', version, 'gbm_purine_denovo_2.pdf', sep = '_'),
       device = "pdf", height = 1, #5.5
       width = 1.5, units = "in")#3.2



# plots for purine simulated MIDs ----

rm(list = ls())
# simulated_data_v9.csv with changed flux bounds paramater_bounds_v3
# simulated_data_v8.csv with Anjali's flux bounds parameter_bounds_v1
version = 'v6'
mid = read.csv(file = paste('../matlab_files_purine/simulated_data_', 
                            version, '.csv', sep = ''), 
               header = T, quote = '')

patient_mids = read.delim(file = '../patient_input_data/patient_MIDs_2405.txt', 
                          header = T, sep = '\t', quote = '')

patient_mids$INFUSION_TIME = patient_mids$INFUSION_TIME/60 # hour

max(patient_mids$INFUSION_TIME) # 3.66 hr
min(patient_mids$INFUSION_TIME) # 2.2 hr

patient_time_points = unique(round(patient_mids$INFUSION_TIME, digits = 1))
patient_time_points

#patient_mids[['psite']] = paste(patient_mids$PATIENT, patient_mids$SITE,sep = '')

patient_mid_rm = subset(patient_mids, 
                        select = -c(R5P0, R5P1, R5P2, R5P3, R5P4, R5P5,
                                    IMP0, GMP0, INO0, GUO0, GDP0, AMP0))

head(patient_mid_rm)

patient_mid_rm$PATIENT = str_sub(patient_mid_rm$PATIENT, start = -1)
patient_mid_rm$SITENO = str_sub(patient_mid_rm$SITE, start = -1)
patient_mid_rm$SITE = str_sub(patient_mid_rm$SITE, start = 1, end = -2)

#patient_mid_rm$SITE = ifelse(str_sub(patient_mid_rm$SITE, start = 1, end = -2) == 'E', 1, 0)

patient_mid_rm$INFUSION_TIME = round(patient_mid_rm$INFUSION_TIME, digits = 1)

patient_mid_long = pivot_longer(patient_mid_rm, cols = 5:ncol(patient_mid_rm)-1)
#patient_mid_rm = patient_mid_rm[, -c(ncol(patient_mid_rm)-1)]
head(patient_mid_long)

patient_mid_long[['PATIENTSITE']] = paste(patient_mid_long$PATIENT, patient_mid_long$SITE, sep = '')

mice_mids = read.delim(file = '../patient_input_data/mice_MMF_gbm_mean_samples.txt',
                       header = T, sep = '\t', quote = '')

mice_mids = subset(mice_mids, select = -c(SER0, SER1, SER2, SER3, 
                                          R5P0, R5P1, R5P2, R5P3, R5P4, R5P5,
                                          IMP0, GMP0, INO0, GUO0, GDP0, AMP0))

mice_mids_long = pivot_longer(mice_mids, cols = 3:ncol(mice_mids))

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

for (i in patient_time_points){
  sim_mid_t = sim_mid_long[sim_mid_long$time == i, ]
  patient_points = patient_mid_long[patient_mid_long$INFUSION_TIME == i, ]
  p = ggplot() + 
    #ggdist::stat_halfeye(adjust = 0.5, justification = -0.3, .width = 0, point_colour = NA, scale = 0.5) +
    geom_violinhalf(data = sim_mid_t, aes(x = name, y = value, fill = Mi), 
                    trim = T, lwd = 0.3, scale = 'width') +
    
    #geom_boxplot(aes(x = fluxNo, y = value),
    #            width = 0.2, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
    geom_half_boxplot(data = sim_mid_t, aes(x = name, y = value),
                      width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
    geom_point(data = patient_points, aes(x = name, y = value, colour = PATIENTSITE),
               shape = 1, size = 0.1) + 
    
    #geom_point(aes(x = 'SER1', y = 0.1), shape = 2, size = 0.1) + 
    
    theme_classic()+
    labs(
      #title = "Raincloud Plot",
      #subtitle = "Showing the Multi-Modal Distributions of Temperatures by Month",
      x    = "Metabolite MID",
      y    = "Simulated MIDs",
      fill = "M+i"
    ) +
    #ggbreak::scale_y_break(c(20, 40)) + # 33 v9, 40 v8
    scale_fill_manual(values = colors) +
    #ylim(-0.5, 30) +
    coord_flip() + 
    theme(legend.key.size = unit(0.05, 'in'), legend.title = element_text(size=7))
  
  #theme(legend.position = 'none')
  
  ggsave(plot = p, filename = paste(version, 'sim_MID', i, '1.pdf', sep = '_'),
         device = "pdf", height = 4, #5.5
         width = 4, units = "in")#3.2
  
}


p = ggplot() + 
  #ggdist::stat_halfeye(adjust = 0.5, justification = -0.3, .width = 0, point_colour = NA, scale = 0.5) +
  geom_violinhalf(data = sim_mid_long, aes(x = name, y = value, fill = Mi), 
                  trim = T, lwd = 0.3, scale = 'width') +
  
  #geom_boxplot(aes(x = fluxNo, y = value),
  #            width = 0.2, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_half_boxplot(data = sim_mid_long, aes(x = name, y = value),
                    width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_point(data = patient_mid_long, aes(x = name, y = value, colour = PATIENTSITE),
             shape = 1, size = 0.1) + 
  #geom_point(aes(x = 'SER1', y = 0.1), shape = 2, size = 0.1) + 
  #geom_point(data = mice_mids_long, aes(x = name, y = value, colour = FileName),
             #shape = 1, size = 0.1) + 
  theme_classic()+
  labs(
    #title = "Raincloud Plot",
    #subtitle = "Showing the Multi-Modal Distributions of Temperatures by Month",
    x    = "Metabolite MID",
    y    = "Simulated MIDs",
    fill = "M+i"
  ) +
  #ggbreak::scale_y_break(c(20, 40)) + # 33 v9, 40 v8
  scale_fill_manual(values = colors) +
  #ylim(-0.5, 20) +
  coord_flip() + 
  theme(legend.key.size = unit(0.05, 'in'), legend.title = element_text(size=7))

#theme(legend.position = 'none')

ggsave(plot = p, 
       filename = paste(version,
                        'sim_MID_all_patients_trim_purine_denovo_gbm.pdf', 
                        sep = '_'),
       device = "pdf", height = 5.5, #5.5
       width = 4, units = "in")#3.2



# tsne plots for purine MID simulation ----
# simulation MIDs
library(Seurat)
library(dplyr)

version = 'v5'
mid = read.csv(file = paste('../matlab_files_purine/glioma_v5/simulated_data_', 
                            version, '.csv', sep = ''), 
               header = T, quote = '')

# simulation parameters including serine and r5p
param_mid = read.csv(file = paste('../matlab_files_purine/glioma_v5/simulation_parameters_', 
                                  version, '.csv', 
                                  sep = ''), header = T, quote = '')

param_mid[['index']] = seq(1, dim(param_mid)[1])

add_mid <- function(mid_df, added_mid, param_df){
  mid_df[[added_mid]] = param_df[match(mid_df$index, param_df$index), added_mid]
  return(mid_df)
}

# combine plasma serine MIDs from parameters to mid dataframe
mid = add_mid(mid, 'ser1', param_mid) 
mid = add_mid(mid, 'ser2', param_mid)  
mid = add_mid(mid, 'ser3', param_mid)  

mid = add_mid(mid, 'r5p1', param_mid) 
mid = add_mid(mid, 'r5p2', param_mid)  
mid = add_mid(mid, 'r5p3', param_mid)  
mid = add_mid(mid, 'r5p4', param_mid) 
mid = add_mid(mid, 'r5p5', param_mid)  

mid[['ser0']] = 1 - mid$ser1 - mid$ser2 - mid$ser3
mid[['r5p0']] = 1 - mid$r5p1 - mid$r5p2 - mid$r5p3 - mid$r5p4 - mid$r5p5

mid = mid[, !(colnames(mid) %in% c('MTHF0', 'MTHF1', 'GLY0', 'GLY1', 'GLY2'))]

mid$time = as.double(mid$time)

mid = subset(mid, time >= 2, select = colnames(mid))

mid_v = mid[, !(colnames(mid) %in% c('index', 'time'))]

# patient MIDs
patient_dir = '../patient_MFA/patient_purine_MFA/patient_mid_mc_gbm_2405'
patient_files = list.files(path = patient_dir)
mid_names = read.delim('../patient_MFA/patient_purine_MFA/mid_name_2405_r.txt', 
                       sep = '\t', header = F)

# combine simulated mids with patient mids
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

pdf(file = 'dimplot_tsne_sim_ps_purine.pdf', height = 10, width = 10)
DimPlot(com_mid, reduction = "tsne", group.by = 'ps_sim', 
        cells = colnames(com_mid)[1501:1051500], cols = my_colors)
DimPlot(com_mid, reduction = "tsne", group.by = 'ps_sim', 
        cells = colnames(com_mid)[1:1500], cols = my_colors)
dev.off()

# imputation ----
version = 'v5'
mid = read.csv(file = paste('../matlab_files_purine/glioma_v5/simulated_data_', 
                            version, '.csv', sep = ''), 
               header = T, quote = '')

# simulation parameters including serine and r5p
param_mid = read.csv(file = paste('../matlab_files_purine/glioma_v5/simulation_parameters_', 
                                  version, '.csv', 
                                  sep = ''), header = T, quote = '')

param_mid[['index']] = seq(1, dim(param_mid)[1])

add_mid <- function(mid_df, added_mid, param_df){
  mid_df[[added_mid]] = param_df[match(mid_df$index, param_df$index), added_mid]
  return(mid_df)
}

# combine plasma serine MIDs from parameters to mid dataframe
mid = add_mid(mid, 'ser1', param_mid) 
mid = add_mid(mid, 'ser2', param_mid)  
mid = add_mid(mid, 'ser3', param_mid)  

mid = add_mid(mid, 'r5p1', param_mid) 
mid = add_mid(mid, 'r5p2', param_mid)  
mid = add_mid(mid, 'r5p3', param_mid)  
mid = add_mid(mid, 'r5p4', param_mid) 
mid = add_mid(mid, 'r5p5', param_mid)  

mid[['ser0']] = 1 - mid$ser1 - mid$ser2 - mid$ser3
mid[['r5p0']] = 1 - mid$r5p1 - mid$r5p2 - mid$r5p3 - mid$r5p4 - mid$r5p5

mid = mid[, !(colnames(mid) %in% c('MTHF0', 'MTHF1', 'GLY0', 'GLY1', 'GLY2',
                                   'ser0', 'r5p0', 'IMP0', 'GMP0', 'GDP0',
                                   'INO0', 'GUO0', 'AMP0'))]

mid$time = as.double(mid$time)

mid = subset(mid, time >= 2, select = colnames(mid))

mid_v = mid[, !(colnames(mid) %in% c('index', 'time'))]

# patient MIDs
patient_dir = '../patient_MFA/patient_purine_MFA/patient_mid_mc_gbm_2405'
patient_files = list.files(path = patient_dir)
mid_names = read.delim('../patient_MFA/patient_purine_MFA/mid_name_2405_r.txt', 
                       sep = '\t', header = F)

# combine simulated mids with patient mids
for (ps in 1:length(patient_files)){
  
  mid_mc_ps = read.delim(paste('../patient_MFA/patient_purine_MFA/patient_mid_mc_gbm_2405/', 
                               patient_files[ps], sep = ''), 
                         header = F, quote = '')
  colnames(mid_mc_ps) = mid_names$V1
  rownames(mid_mc_ps) = paste(rep(paste('P', ps, sep = '')), 
                              seq(1, dim(mid_mc_ps)[1]), sep = '_')
  
  mid_mc_ps = mid_mc_ps[, !colnames(mid_mc_ps) %in% c('ser0', 'r5p0', 
                                                      'IMP0', 'GMP0', 'GDP0',
                                                      'INO0', 'GUO0', 'AMP0')]
  
  mid_v = rbind(mid_mc_ps, mid_v)
}

mid_v = mid_v * 100
mid_v = data.matrix(mid_v)
mid_v = t(mid_v)

# gbm 5e4
sim_mid = mid_v[, 1501:1051500]
patient_mid = mid_v[, 1:1500]

# gbm 1e4
sim_mid = mid_v[, 1501:211500]
patient_mid = mid_v[, 1:1500]
#cortex
sim_mid = mid_v[, 801:210800]
patient_mid = mid_v[, 1:800]

sim_mid = CreateSeuratObject(counts = sim_mid, project = 'sim_data_purine_v5_1e4')
patient_mid = CreateSeuratObject(counts = patient_mid, project = 'patient_mid_mc')

saveRDS(sim_mid, 'sim_gbm_purine_mid_v5_1e4.rds')
saveRDS(patient_mid, 'sim_patient_purine_mid_mc.rds')
  

sim_mid = SCTransform(sim_mid) %>% RunPCA()
sim_mid = RunUMAP(sim_mid, dims = 1:30)

patient_mid = SCTransform(patient_mid) %>% RunPCA()
anchors = FindTransferAnchors(reference = sim_mid, 
                              query = patient_mid, 
                              normalization.method = "SCT")

predictions.assay = TransferData(anchorset = anchors, 
                                 refdata = GetAssayData(sim_mid[['RNA']]),
                                 prediction.assay = T,
                                 weight.reduction = patient_mid[["pca"]],
                                 dims = 1:30)

pred_assay = as.data.frame(t(as.matrix(predictions.assay@data)))

write.csv(pred_assay, row.names = T, file = 'patient_purine_gbm_2405_5e4_rm0_100.csv')

for (i in 1:15){
  f = (i - 1) * 100 + 1
  e = i * 100
  write.csv(pred_assay[f:e, ], row.names = T,
            file = paste(substr(rownames(pred_assay)[i * 100], 1, 3),
                         'patient_purine_gbm_2405_5e4_rm0_100.csv', sep = '_'))
}

saveRDS(sim_mid, 'sim_cortex_purine_mid_v5_1e4_rm0_100_sct.rds')
saveRDS(patient_mid, 'sim_patient_purine_mid_mc_sct.rds')

sim_mid = readRDS('sim_gbm_purine_mid_v5_1e4_rm0_100_sct.rds')

# mice MIDs
mice_dir = '../patient_MFA/patient_purine_MFA/mice_mid_mc_2406'
mice_dir = '../patient_MFA/patient_purine_MFA/TRP_GBM38_mid_mc/'
mice_files = list.files(path = mice_dir)
mid_names = read.delim('../patient_MFA/patient_purine_MFA/mid_name_2405_r.txt', 
                       sep = '\t', header = F)

mice_mid = c()
for (ps in 1:length(mice_files)){
  
  mid_mc_ps = read.delim(paste(mice_dir, 
                               mice_files[ps], sep = ''), 
                         header = F, quote = '')
  colnames(mid_mc_ps) = mid_names$V1
  rownames(mid_mc_ps) = paste(rep(paste('P', ps, sep = '')), 
                              seq(1, dim(mid_mc_ps)[1]), sep = '_')
  
  mid_mc_ps = mid_mc_ps[, !colnames(mid_mc_ps) %in% c('ser0', 'r5p0', 
                                                      'IMP0', 'GMP0', 'GDP0',
                                                      'INO0', 'GUO0', 'AMP0')]
  
  mice_mid = c(mice_mid, list(mid_mc_ps))
}
names(mice_mid) = substr(mice_files, 1, 20)
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
write.csv(pred_mice_assay, row.names = T, file = 'mice_purine_gbm_2405_5e4_rm0_100.csv')
write.csv(pred_mice_assay, row.names = T, file = 'TRP_GBM38_1e4_rm0_100_2.csv')

for (i in 1:4){ #20
  f = (i - 1) * 100 + 1
  e = i * 100
  write.csv(pred_mice_assay[f:e, ], row.names = T,
            file = paste(names(mice_mid)[i],
                         'mice_purine_gbm_2406_1e4_rm0_100.csv', sep = '_'))
}

# mice MIDs - cortex
mice_dir = '../patient_MFA/patient_purine_MFA/mice_mid_mc_cortex_2405'
mice_files = list.files(path = mice_dir)
mid_names = read.delim('../patient_MFA/patient_purine_MFA/mid_name_2405_r.txt', 
                       sep = '\t', header = F)

mice_mid = c()
for (ps in 1:length(mice_files)){
  
  mid_mc_ps = read.delim(paste('../patient_MFA/patient_purine_MFA/mice_mid_mc_cortex_2405/', 
                               mice_files[ps], sep = ''), 
                         header = F, quote = '')
  colnames(mid_mc_ps) = mid_names$V1
  rownames(mid_mc_ps) = paste(rep(paste('P', ps, sep = '')), 
                              seq(1, dim(mid_mc_ps)[1]), sep = '_')
  
  mice_mid = c(mice_mid, list(mid_mc_ps))
}
names(mice_mid) = substr(mice_files, 1, 10)
mice_mid_v = dplyr::bind_rows(mice_mid)

mice_mid_v = data.matrix(mice_mid_v)
mice_mid_v = t(mice_mid_v)

mice_mid_v = CreateSeuratObject(counts = mice_mid_v, project = 'mice_mid_mc_cortex')

saveRDS(mice_mid_v, 'sim_mice_mid_mc_cortex_2405.rds')

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
write.csv(pred_mice_assay, row.names = T, file = 'mice_purine_cortex_2405.csv')

for (i in 1:20){
  f = (i - 1) * 100 + 1
  e = i * 100
  write.csv(pred_mice_assay[f:e, ], row.names = T,
            file = paste(names(mice_mid)[i],
                         'mice_purine_cortex_2405.csv', sep = '_'))
}

# plot MMF mice gbm de novo GMP prediction ----
#mice_pred_dir = '../../ML_serine_model/metabolic_CNN/mice_mmf_2406_1e4_rm0_100_240517_1138'
#mice_pred_dir = '../../ML_serine_model/mice_mmf_cortex_240521_1854_imputed_1e4_rm0_100'
mice_pred_dir = '../../ML_serine_model/mice_mmf_glioma_240517_1138_imputed_1e4_rm0_100'
#mice_pred_dir = '../../ML_serine_model/mice_mmf_cortex_240521_1854_imputed_2405'
mice_pred_files = list.files(path = mice_pred_dir)

list_mice_pred = c()
for (ps in 1:length(mice_pred_files)){
  mice_pred = read.delim(file = paste(mice_pred_dir, mice_pred_files[ps], 
                                      sep = '/'), 
                         sep = '\t', header = T)
  
  list_mice_pred = c(list_mice_pred, list(mice_pred$denovo))
}

names(list_mice_pred) = substr(mice_pred_files, 42, 50)  #gbm
#names(list_mice_pred) = substr(mice_pred_files, 30, 38) #cortex

mice_pred_df = dplyr::bind_rows(list_mice_pred)
drops <- c("high2_773", "high2_774", "high2_780",
           "high4_775", "high4_778", "high4_779",
           "low2_766G", "low2_767G", "low2_782G",
           "low4_769G", "low4_771G", "low4_772G")

drops <- c("high1_764", "high1_765", "high1_781", 
           "high2_773", "high2_774", "high2_780",
           "high4_775", "high4_778", "high4_779")

mice_pred_df = mice_pred_df[, !(names(mice_pred_df) %in% drops)]


mice_pred_long = tidyr::pivot_longer(mice_pred_df, cols = 1:ncol(mice_pred_df))

treatment_sample = stringr::str_split_fixed(mice_pred_long$name, "_", 2)

treatment_sample = as.data.frame(treatment_sample)
mice_pred_long[['treatment']] = treatment_sample$V1
mice_pred_long[['sample']] = treatment_sample$V2

ReplicateAverages = mice_pred_long %>% 
  group_by(treatment, sample) %>% 
  summarise_each(list(mean)) 

TreatmentAverages = mice_pred_long %>% 
  group_by(treatment) %>% 
  summarise_each(list(mean)) 

ReplicateAverages.summary <- ReplicateAverages %>%
  group_by(treatment) %>%
  summarise(sd = sd(value, na.rm = T),
            value = mean(value))

ReplicateAverages.summary

p = ggplot(ReplicateAverages, aes(treatment, value)) +
  geom_bar(stat = "identity", data = ReplicateAverages.summary,
           fill = NA, color = 'black') +
  geom_jitter(position = position_jitter(0.2),
              color = 'black') +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),
                data = ReplicateAverages.summary, width = 0.2)

p
#https://rupress.org/jcb/article/219/6/e202001064/151717/SuperPlots-Communicating-reproducibility-and
mice_pred_long_2 = mice_pred_long[, -c(1)]
mice_pred_long_3 = mice_pred_long_2[, -c(3)]

mice_pred_long_3$treatment = factor(mice_pred_long_3$treatment, 
                                    levels = c('ctrl', 'low1', 'low2', 'low4'))

mice_pred_long_2$treatment = factor(mice_pred_long_2$treatment, 
                                    levels = c('ctrl', 'low1', 'high1'))

ReplicateAverages = mice_pred_long_2 %>% 
  group_by(treatment, sample) %>% 
  summarise_each(list(mean)) 

TreatmentAverages = mice_pred_long_3 %>% 
  group_by(treatment) %>% 
  summarise_each(list(mean)) 

#ReplicateAverages[['sample']] = NULL

p = ggplot(mice_pred_long_2, aes(x=treatment,y=value,color=factor(sample))) + 
  #geom_jitter() +
  geom_beeswarm(corral = "wrap", corral.width = 0.9, cex = 1, size = 0.5) + 
  #scale_colour_brewer(palette = "Set3") + 
  scale_color_manual(values = c('#90AD1C', #CTRL
                                '#C075A6', #LOW1
                                '#AA0DFE', #LOW1
                                '#1C8356', #CTRL
                                '#DEA0FD', #LOW1
                                '#325A9B', #LOW2
                                '#2ED9FF', #LOW2
                                '#C4451C', #LOW4
                                '#FEAF16', #LOW4 
                                '#822E1C', #LOW4
                                '#3283FE' #LOW2
                                )) +
  stat_compare_means(data = ReplicateAverages,
                     comparisons = list(c('ctrl', 'low1'),
                                        c('ctrl', 'low2'),
                                        c('ctrl', 'low4'),
                                        c('low1', 'low2'),
                                        c('low1', 'low4'),
                                        c('low2', 'low4')),
                     #method = 'anova',
                     label.y = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.95),
                     size = 2) + 
  geom_beeswarm(data = ReplicateAverages, size = 3, shape = 15) + 
  geom_beeswarm(data = TreatmentAverages, size = 3, shape = 15, color = '#474747') +
  theme_classic()
  #theme(legend.position="none")

ggsave(plot = p, filename = 'mice_mmf_glioma_240517_1138_imputed_5e4_rm0_100_doses.pdf',
       device = "pdf", height = 5, #5.5
       width = 6, units = "in")#3.2


p = ggplot(mice_pred_long_2, aes(x=treatment,y=value,color=factor(sample))) + 
  #geom_jitter() +
  geom_beeswarm(corral = "wrap", corral.width = 0.9, cex = 1, size = 0.5) + 
  #scale_colour_brewer(palette = "Set3") + 
  scale_color_manual(values = c('#90AD1C', #CTRL
                                '#C075A6', #LOW1
                                '#AA0DFE', #LOW1
                                '#1C8356', #CTRL
                                '#DEA0FD', #LOW1
                                '#C4451C', #high1
                                '#FEAF16', #high1
                                '#822E1C' #high1
                                )) +
  stat_compare_means(#data = ReplicateAverages,
    comparisons = list(c('ctrl', 'low1'),
                       c('ctrl', 'high1'),
                       c('low1', 'high1')),
    #method = 'anova',
    label.y = c(0.05, 0.1, 0.95),
    size = 2) + 
  geom_beeswarm(data = ReplicateAverages, size = 3, shape = 15) + 
  geom_beeswarm(data = TreatmentAverages, size = 3, shape = 15, color = '#474747') +
  theme_classic()
#theme(legend.position="none")

ggsave(plot = p, filename = 'mice_mmf_glioma_240517_1138_imputed_5e4_rm0_100_conc.pdf',
       device = "pdf", height = 5, #5.5
       width = 6, units = "in")#3.2



# 2406
names(list_mice_pred) = c('ctrl', 'low1', 'low2', 'low4')

mice_pred_df = dplyr::bind_rows(list_mice_pred)

mice_pred_long = tidyr::pivot_longer(mice_pred_df, cols = 1:ncol(mice_pred_df))

df.summary <- mice_pred_long %>%
  group_by(name) %>%
  summarise(sd = sd(value, na.rm = T),
            value = mean(value))
df.summary

p = ggplot(mice_pred_long, aes(name, value)) +
  geom_bar(stat = "identity", data = df.summary,
           fill = NA, color = 'black') +
  geom_jitter(position = position_jitter(0.2),
              color = 'black') +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),
                data = df.summary, width = 0.2)

p

p = ggplot(mice_pred_long, aes(name, value)) +
  geom_violin(trim = F, scale = 'width')


a1 = aov(data = ReplicateAverages, formula = value ~ treatment)
pairwise.t.test(ReplicateAverages$value, ReplicateAverages$treatment, 
                p.adj = "bonf")
summary(a1)
TukeyHSD(a1)

