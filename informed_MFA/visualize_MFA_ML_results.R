library(RColorBrewer)
library(see)
library(ggplot2)
library(dplyr)
library(gghalves)
library(tidyr)
library(readxl)

# plot MFA-estimated cortex fluxes relative to CNN-predicted glioma fluxes ----
# Figure 3K-L

# MFA-estimated cortex fluxes relative to CNN predicted glioma fluxes
sample_dir = list.files(path = './output_files/', 
                        pattern = 'Patient', all.files = T, full.names = T)

sample_names = strsplit(sample_dir, '/')
sample_names = lapply(sample_names, function (x) x[3])
sample_names = unlist(sample_names)

sample_names = strsplit(sample_names, '_')
sample_names = lapply(sample_names, function (x) x[1])
sample_names = unlist(sample_names)

denovo_cortex = c()
plasma_cortex = c()

for (ps in 1:length(sample_names)){

  flux = readxl::read_excel(path = paste(sample_dir[ps], 'flux_results_ML.xlsx', 
                                         sep = '/'), sheet = 'flux') 
  denovo_c = flux[flux$Reaction == 'PGa == SERa', c(2:dim(flux)[2])]
  plasma_c = flux[flux$Reaction == 'SERx == SERa', c(2:dim(flux)[2])]
                            
  denovo_cortex = c(denovo_cortex, list(t(denovo_c)))
  plasma_cortex = c(plasma_cortex, list(t(plasma_c)))
} 

denovo_cortex = bind_cols(denovo_cortex)
plasma_cortex = bind_cols(plasma_cortex)

colnames(denovo_cortex) = sample_names
colnames(plasma_cortex) = sample_names

# CNN-predicted glioma fluxes
ML_glioma_plasma_dir = list.files(path = '../metabolic_CNN/serine_plasma_glioma_config_20240408_2102/patients_pred_serine_plasma_glioma', 
                                  pattern = 'P', all.files = T, full.names = T)
ML_glioma_denovo_dir = list.files(path = '../metabolic_CNN/serine_denovo_glioma_config_20240408_1327/patients_pred_serine_denovo_glioma', 
                                  pattern = 'P', all.files = T, full.names = T)

denovo_glioma = c()
plasma_glioma = c()

for (ps in 1:length(sample_names)){
    
  ML_ratio_glioma_plasma = read.delim(file = ML_glioma_plasma_dir[ps], 
                                      quote = '', header = T)
  ML_ratio_glioma_denovo = read.delim(file = ML_glioma_denovo_dir[ps], 
                                      quote = '', header = T)
  
  denovo_glioma = c(denovo_glioma, as.data.frame(ML_ratio_glioma_denovo$denovo))
  plasma_glioma = c(plasma_glioma, as.data.frame(ML_ratio_glioma_plasma$plasma))
  
}

names(denovo_glioma) = sample_names
names(plasma_glioma) = sample_names

denovo_glioma = bind_cols(denovo_glioma)
plasma_glioma = bind_cols(plasma_glioma)


combined_MFA_ML_denovo = rbind(denovo_cortex, denovo_glioma)
combined_MFA_ML_plasma = rbind(plasma_cortex, plasma_glioma)

combined_MFA_ML_denovo[['cortex_glioma']] = c(rep('cortex', 100), rep('glioma', 100))
combined_MFA_ML_plasma[['cortex_glioma']] = c(rep('cortex', 100), rep('glioma', 100))

denovo_long = tidyr::gather(combined_MFA_ML_denovo, ps_id, value, Patient1E:Patient8N, 
                     factor_key = T)
plasma_long = tidyr::gather(combined_MFA_ML_plasma, ps_id, value, Patient1E:Patient8N, 
                     factor_key = T)

denovo_long[['patient']] = substring(as.character(denovo_long$ps_id), 1, 8)
plasma_long[['patient']] = substring(as.character(plasma_long$ps_id), 1, 8)


colours <- c("#50C878", "#ff8007") 

# replace denovo_long with plasma_long to create Figure 3L
p = ggplot(denovo_long, aes(x = patient, y = value, fill = cortex_glioma)) +
  introdataviz::geom_split_violin(alpha = 0.6, scale = 'width') +
  geom_boxplot(outlier.shape = NA, width = 0.3,
               show.legend = FALSE) +
  scale_x_discrete(name = "Patient") +
  scale_y_continuous(name = "Glucose-derived serine synthesis flux") +
  scale_fill_manual(values = colours)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(filename = 'vlnplot_cortex_vs_glioma_denovo_MFA_ML.pdf', plot = p, 
       device = 'pdf', height = 4, width = 5, units = 'in')

