library(RColorBrewer)
library(see)
library(ggplot2)
library(dplyr)
library(gghalves)
library(tidyr)

rm(list = ls())
version = 'v29'
flux = read.delim(file = paste('../matlab_files_v2/flux_', version, '.txt', sep = ''),
                  header = F, sep = '\t', quote = '')

bounds = read.delim(file = paste('../matlab_files_v2/flux_bounds_', version, '.txt', sep = ''),
                    header = F, sep = '\t', quote = '')

rxn_ids = read.delim(file = paste('../matlab_files_v2/rxn_ids_', version, '.txt', sep = ''),
                     header = F, sep = '\t', quote = '')

no_rxns = nrow(rxn_ids)
flux = t(flux)
head(rownames(flux))
dim(flux)
colnames(flux) = paste(rep('v', no_rxns), seq(1, no_rxns), sep = '')

bounds[['reactions']] = rxn_ids$V1
bounds[['flux_no']] = colnames(flux)
colnames(bounds) = c('lb', 'ub', 'reactions', 'flux_no')

flux_names = paste(rep('v', no_rxns), seq(1, no_rxns), sep = '')
flux = as.data.frame(flux)
flux[['sim']] = factor(rownames(flux))

flux_long = gather(flux, fluxNo, value, v1:v55, factor_key = T)

net_rxn = c('v1', 'v3', 'v5', 'v7', 'v9', 'v11', 'v13', 'v15', 'v17')
net_rxn_ratio = c('v2', 'v4', 'v6', 'v8', 'v10', 'v12', 'v14', 'v16', 'v18')
irrev_rxn = paste(rep('v', no_rxns - 19 + 1), seq(19, no_rxns), sep = '')
exchange_rxn = c('v34', 'v35', 'v36', 'v37', 'v38', 'v39', 'v40', 
                 'v46', 'v49', 'v54', 'v55')
#exchange_rxn = c('v32', 'v33', 'v34', 'v35', 'v36', 'v37', 'v38', 
                 #'v44', 'v47', 'v52', 'v53')
irrev_rxn = irrev_rxn[!irrev_rxn %in% exchange_rxn]

hex_codes = c()
for (x in flux_names){
  flux_long[flux_long$fluxNo == x, 'lb'] = bounds[bounds$flux_no == x, 'lb']  
  flux_long[flux_long$fluxNo == x, 'ub'] = bounds[bounds$flux_no == x, 'ub']  
  
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
                  trim = T, lwd = 0.3, scale = 'width') +

  geom_half_boxplot(aes(x = fluxNo, y = value),
                    width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_point(aes(x = fluxNo, y = lb), shape = 6, size = 0.1) + 
  geom_point(aes(x = fluxNo, y = ub), shape = 2, size = 0.1) + 
  
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

ggsave(plot = p, filename = paste('flux_distribution', version, '.jpg', sep = '_'),
       device = "jpeg", height = 6, #5.5
       width = 3.2, units = "in", dpi = 1000)#3.2

p = flux_long_selected_2 %>%
  ggplot() + 
  geom_violinhalf(aes(x = fluxNo, y = value, fill = fluxNo), 
                  trim = T, lwd = 0.3, scale = 'width') +
  
  geom_half_boxplot(aes(x = fluxNo, y = value),
                    width = 0.6, outlier.color = NA, alpha = 0.8, lwd = 0.3) + 
  geom_point(aes(x = fluxNo, y = lb), shape = 6, size = 0.1) + 
  geom_point(aes(x = fluxNo, y = ub), shape = 2, size = 0.1) + 
  
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

ggsave(plot = p, filename = paste('flux_distribution', version, '2.jpg', sep = '_'),
       device = "jpeg", height = 2, #5.5
       width = 2.5, units = "in", dpi = 1000)#3.2

###############################################################################
# plot MFA results vs ML results 

# read MFA results
sample_dir = list.files(path = '../patient_MFA/patient_serine_MFA/Output files/serine_MFA_second_try_scflux', 
                        pattern = 'Patient', all.files = T, full.names = T)

sample_names = strsplit(sample_dir, '/')
sample_names = lapply(sample_names, function (x) x[6])
sample_names = unlist(sample_names)

sample_names = strsplit(sample_names, '_')
sample_names = lapply(sample_names, function (x) x[1])
sample_names = unlist(sample_names)

# read ML results for cortex
ML_cortex_dir = list.files(path = '../../ML_serine_model/pred_patient_cortex_20240403_1615', 
                        pattern = 'pred_patient_denovo_cortex', all.files = T, full.names = T)


# MFA cortex
ratio4_appended = c()
for (ps in 1:length(sample_names)){
  #ratio_bounds = read.delim(file = paste(sample_dir[ps], 'flux_ratio.txt', sep = '/'), 
                          #quote = '', header = 1, row.names = 1)

  ratio_cortex = read.delim(file = paste(sample_dir[ps], 'ratio4_cortex.txt', 
                                         sep = '/'), 
                            quote = '', header = F)
  ratio4_appended = c(ratio4_appended, ratio_cortex)
}  

# ML cortex
ML_ratio4_appended = c()
ML_cortex_plasma_appended = c()
for (ps in 1:length(sample_names)){
  #ratio_bounds = read.delim(file = paste(sample_dir[ps], 'flux_ratio.txt', sep = '/'), 
  #quote = '', header = 1, row.names = 1)
  
  ML_ratio_cortex = read.delim(file = ML_cortex_dir[ps],
                            quote = '', header = T)
  ML_ratio4_appended = c(ML_ratio4_appended, as.data.frame(ML_ratio_cortex$denovo))
  ML_cortex_plasma_appended = c(ML_cortex_plasma_appended, 
                                as.data.frame(1 - ML_ratio_cortex$denovo))
}  

names(ratio4_appended) = sample_names
names(ML_ratio4_appended) = sample_names
names(ML_cortex_plasma_appended) = sample_names

MFA_ratios = bind_cols(ratio4_appended)
ML_ratios = bind_cols(ML_ratio4_appended)

combined_MFA_ML = rbind(MFA_ratios, ML_ratios)
combined_MFA_ML[['MFA_ML']] = c(rep('MFA', 100), rep('ML', 100))

ratio_long = gather(combined_MFA_ML, ps_id, value, Patient1E:Patient8N, 
                    factor_key = T)

ratio_long[['patient']] = substring(as.character(ratio_long$ps_id), 1, 8)


colours <- c("#ff8007", "#05badd")

p = ggplot(ratio_long, aes(x = patient, y = value, fill = MFA_ML)) +
  introdataviz::geom_split_violin(alpha = 0.6, scale = 'width') +
  geom_boxplot(outlier.shape = NA, width = 0.3,
               show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = 'red') +
  scale_x_discrete(name = "Patient") +
  scale_y_continuous(name = "Glucose-derived serine ratio (cortex)") +
  scale_fill_manual(values = colours)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(filename = 'vlnplot_patient_cortex_MFA_ML.pdf', plot = p, 
       device = 'pdf', height = 4, width = 6, units = 'in')

# MFA glioma
ratio_appended_1 = c()
ratio_appended_2 = c()
ratio_appended_3 = c()
for (ps in 1:length(sample_names)){
  #ratio_bounds = read.delim(file = paste(sample_dir[ps], 'flux_ratio.txt', sep = '/'), 
  #quote = '', header = 1, row.names = 1)
  
  ratio_glioma = read.delim(file = paste(sample_dir[ps], 'xmc_glioma.txt', 
                                         sep = '/'), 
                            quote = '', header = F)
  ratio_appended_1 = c(ratio_appended_1, as.data.frame(ratio_glioma$V1))
  ratio_appended_2 = c(ratio_appended_2, as.data.frame(ratio_glioma$V2))
  ratio_appended_3 = c(ratio_appended_3, as.data.frame(ratio_glioma$V3))
}

# ML glioma
ML_glioma_plasma_dir = list.files(path = '../../ML_serine_model/pred_patient_glioma_plasma_20240408_2102', 
                                   pattern = 'pred_patient_glioma', all.files = T, full.names = T)
ML_glioma_denovo_dir = list.files(path = '../../ML_serine_model/pred_patient_glioma_denovo_20240408_1327', 
  pattern = 'pred_patient_glioma', all.files = T, full.names = T)

ML_ratio_appended_1 = c()
ML_ratio_appended_2 = c()
ML_ratio_appended_3 = c()
for (ps in 1:length(sample_names)){
  #ratio_bounds = read.delim(file = paste(sample_dir[ps], 'flux_ratio.txt', sep = '/'), 
  #quote = '', header = 1, row.names = 1)
  
  ML_ratio_glioma_plasma = read.delim(file = ML_glioma_plasma_dir[ps], 
                               quote = '', header = T)
  ML_ratio_glioma_denovo = read.delim(file = ML_glioma_denovo_dir[ps], 
                                      quote = '', header = T)
  
  ML_ratio_appended_1 = c(ML_ratio_appended_1, as.data.frame(ML_ratio_glioma_denovo$denovo))
  ML_ratio_appended_2 = c(ML_ratio_appended_2, as.data.frame(ML_ratio_glioma_plasma$plasma))
  ML_ratio_appended_3 = c(ML_ratio_appended_3, 
                          as.data.frame(1 - ML_ratio_glioma_denovo$denovo - ML_ratio_glioma_plasma$plasma))
}



                                  
ML_glioma_dir = list.files(path = '../../ML_serine_model/pred_patient_glioma_tme_20240409_1120', 
                           pattern = 'pred_patient_glioma', all.files = T, full.names = T)
#ML_ratio_appended_1 = c()
#ML_ratio_appended_2 = c()
ML_ratio_appended_3 = c()
for (ps in 1:length(sample_names)){
  #ratio_bounds = read.delim(file = paste(sample_dir[ps], 'flux_ratio.txt', sep = '/'), 
  #quote = '', header = 1, row.names = 1)
  
  ML_ratio_glioma = read.delim(file = ML_glioma_dir[ps], 
                            quote = '', header = T)
  #ML_ratio_appended_1 = c(ML_ratio_appended_1, as.data.frame(ML_ratio_glioma$denovo))
  #ML_ratio_appended_2 = c(ML_ratio_appended_2, as.data.frame(ML_ratio_glioma$plasma))
  ML_ratio_appended_3 = c(ML_ratio_appended_3, as.data.frame(ML_ratio_glioma$tme))
}


plot_glioma_MFA_ML = function(ratio_appended, ML_ratio_appended, filename){
  names(ratio_appended) = sample_names  
  names(ML_ratio_appended) = sample_names 
  
  MFA_ratios = bind_cols(ratio_appended)
  ML_ratios = bind_cols(ML_ratio_appended)
  
  combined_MFA_ML = rbind(MFA_ratios, ML_ratios)
  combined_MFA_ML[['MFA_ML']] = c(rep('MFA', 100), rep('ML', 100))
  
  ratio_long = tidyr::gather(combined_MFA_ML, ps_id, value, Patient1E:Patient8N, 
                      factor_key = T)
  
  ratio_long[['patient']] = substring(as.character(ratio_long$ps_id), 1, 8)
  
  
  colours <- c("#ff8007", "#05badd")
  
  p = ggplot(ratio_long, aes(x = ps_id, y = value, fill = MFA_ML)) +
    introdataviz::geom_split_violin(alpha = 0.6, scale = 'width') +
    geom_boxplot(outlier.shape = NA, width = 0.2,
                 show.legend = FALSE) +
    scale_x_discrete(name = "Patient") +
    scale_y_continuous(name = filename) +
    scale_fill_manual(values = colours)+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  
  ggsave(filename = paste('vlnplot_patient_glioma_MFA_ML_combined_1327_2102', filename, '.pdf', sep = ''), 
         plot = p, 
         device = 'pdf', height = 4, width = 8, units = 'in')
}

plot_glioma_MFA_ML(ratio_appended_1, ML_ratio_appended_1, 'Glucose-derived serine')
plot_glioma_MFA_ML(ratio_appended_2, ML_ratio_appended_2, 'Plasma serine')
plot_glioma_MFA_ML(ratio_appended_3, ML_ratio_appended_3, 'TME serine')


names(ML_ratio_appended_1) = sample_names
names(ML_ratio_appended_2) = sample_names
names(ML_ratio_appended_3) = sample_names

ML_ratio_appended_1 = as.data.frame(ML_ratio_appended_1)
ML_ratio_appended_2 = as.data.frame(ML_ratio_appended_2)
ML_ratio_appended_3 = as.data.frame(ML_ratio_appended_3)

write.table(ML_ratio_appended_1, file = 'ML_denovo_combined_1327_2102.txt', 
            sep = "\t", quote = F,
            row.names = F, col.names = T)
write.table(ML_ratio_appended_2, file = 'ML_plasma_combined_1327_2102.txt', 
            sep = "\t", quote = F,
            row.names = F, col.names = T)
write.table(ML_ratio_appended_3, file = 'ML_tme_combined_1327_2102.txt', 
            sep = "\t", quote = F,
            row.names = F, col.names = T)
write.table(ML_ratio4_appended, file = 'ML_denovo_cortex_1615.txt', 
            sep = "\t", quote = F,
            row.names = F, col.names = T)

write.table(ML_cortex_plasma_appended, file = 'ML_plasma_cortex_1615.txt', 
            sep = "\t", quote = F,
            row.names = F, col.names = T)

###############################################################################
# plot MFA ML
library(readxl)

# read MFA results
sample_dir = list.files(path = '../patient_MFA/patient_serine_MFA/Output files/serine_MFA_ML_separated_models', 
                        pattern = 'Patient', all.files = T, full.names = T)

sample_names = strsplit(sample_dir, '/')
sample_names = lapply(sample_names, function (x) x[6])
sample_names = unlist(sample_names)

sample_names = strsplit(sample_names, '_')
sample_names = lapply(sample_names, function (x) x[1])
sample_names = unlist(sample_names)

dc = c()
pc = c()
for (ps in 1:length(sample_names)){
  #ratio_bounds = read.delim(file = paste(sample_dir[ps], 'flux_ratio.txt', sep = '/'), 
  #quote = '', header = 1, row.names = 1)
  
  flux = read_excel(path = paste(sample_dir[ps], 'flux_results_ML.xlsx', 
                                         sep = '/'), sheet = 'flux') 
  denovo_c = flux[flux$Reaction == 'PGa == SERa', c(2:dim(flux)[2])]
  plasma_c = flux[flux$Reaction == 'SERx == SERa', c(2:dim(flux)[2])]
                            
  dc = c(dc, list(t(denovo_c)))
  pc = c(pc, list(t(plasma_c)))
} 

dc = bind_cols(dc)
pc = bind_cols(pc)

colnames(dc) = sample_names
colnames(pc) = sample_names

# ML glioma
ML_glioma_plasma_dir = list.files(path = '../../ML_serine_model/pred_patient_glioma_plasma_20240408_2102', 
                                  pattern = 'pred_patient_glioma', all.files = T, full.names = T)
ML_glioma_denovo_dir = list.files(path = '../../ML_serine_model/pred_patient_glioma_denovo_20240408_1327', 
                                  pattern = 'pred_patient_glioma', all.files = T, full.names = T)

ML_ratio_appended_1 = c()
ML_ratio_appended_2 = c()

for (ps in 1:length(sample_names)){
  #ratio_bounds = read.delim(file = paste(sample_dir[ps], 'flux_ratio.txt', sep = '/'), 
  #quote = '', header = 1, row.names = 1)
  
  ML_ratio_glioma_plasma = read.delim(file = ML_glioma_plasma_dir[ps], 
                                      quote = '', header = T)
  ML_ratio_glioma_denovo = read.delim(file = ML_glioma_denovo_dir[ps], 
                                      quote = '', header = T)
  
  ML_ratio_appended_1 = c(ML_ratio_appended_1, as.data.frame(ML_ratio_glioma_denovo$denovo))
  ML_ratio_appended_2 = c(ML_ratio_appended_2, as.data.frame(ML_ratio_glioma_plasma$plasma))
  
}

names(ML_ratio_appended_1) = sample_names
names(ML_ratio_appended_2) = sample_names

ML_ratio_appended_1 = bind_cols(ML_ratio_appended_1)
ML_ratio_appended_2 = bind_cols(ML_ratio_appended_2)


combined_MFA_ML_denovo = rbind(dc, ML_ratio_appended_1)
combined_MFA_ML_plasma = rbind(pc, ML_ratio_appended_2)

combined_MFA_ML_denovo[['MFA_ML']] = c(rep('MFA', 100), rep('ML', 100))
combined_MFA_ML_plasma[['MFA_ML']] = c(rep('MFA', 100), rep('ML', 100))

denovo_long = tidyr::gather(combined_MFA_ML_denovo, ps_id, value, Patient1E:Patient8N, 
                     factor_key = T)
plasma_long = tidyr::gather(combined_MFA_ML_plasma, ps_id, value, Patient1E:Patient8N, 
                     factor_key = T)

denovo_long[['patient']] = substring(as.character(denovo_long$ps_id), 1, 8)
plasma_long[['patient']] = substring(as.character(plasma_long$ps_id), 1, 8)


colours <- c("#50C878", "#ff8007") #"#05badd"

p = ggplot(denovo_long, aes(x = patient, y = value, fill = MFA_ML)) +
  introdataviz::geom_split_violin(alpha = 0.6, scale = 'width') +
  geom_boxplot(outlier.shape = NA, width = 0.3,
               show.legend = FALSE) +
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = 'red') +
  scale_x_discrete(name = "Patient") +
  scale_y_continuous(name = "Glucose-derived serine synthesis flux") +
  scale_fill_manual(values = colours)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(filename = 'vlnplot_cortex_vs_glioma_denovo_MFA_ML.pdf', plot = p, 
       device = 'pdf', height = 4, width = 5, units = 'in')


###############################################################################
# plot M+3 serine to M+3 PG and M+1 serine to M+1 plasma serine
ML_pred_dir = "C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/data_simulation/plot_scripts/separated_models_2404/"

glioma_denovo = read.delim(file = paste(ML_pred_dir, 'ML_denovo_combined_1327_2102.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')

glioma_plasma = read.delim(file = paste(ML_pred_dir, 'ML_plasma_combined_1327_2102.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
glioma_tme = read.delim(file = paste(ML_pred_dir, 'ML_tme_combined_1327_2102.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
cortex_denovo = read.delim(file = paste(ML_pred_dir, 'ML_denovo_cortex_1615.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
cortex_plasma = read.delim(file = paste(ML_pred_dir, 'ML_plasma_cortex_1615.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')

patient_site_names = colnames(glioma_denovo)

patient_folders = list.dirs(path = "C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/data_simulation/patient_MFA/patient_serine_MFA/Output files/serine_MFA_second_try_scflux",
                            recursive = F)
mid_names = read.delim(file = "C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/data_simulation/patient_MFA/patient_serine_MFA/mid_name.txt", 
                       header = F, sep = '\t')

for (ps in 1:ncol(glioma_denovo)){
patient_mids = read.delim(file = paste(patient_folders[ps], '/mid_mc.txt', sep = ''), 
           header = F, sep = '\t')


colnames(patient_mids) = mid_names$V1

glioma_denovo[[paste(patient_site_names[ps], 'M3', sep = '_')]] = patient_mids$SER3 / patient_mids$PG3
glioma_plasma[[paste(patient_site_names[ps], 'M1', sep = '_')]] = patient_mids$SER1 / patient_mids$SERp1

glioma_tme[[paste(patient_site_names[ps], 'M1', sep = '_')]] = patient_mids$SER1 / patient_mids$SERc1
glioma_tme[[paste(patient_site_names[ps], 'M2', sep = '_')]] = patient_mids$SER2 / patient_mids$SERc2
glioma_tme[[paste(patient_site_names[ps], 'M3', sep = '_')]] = patient_mids$SER3 / patient_mids$SERc3

cortex_denovo[[paste(patient_site_names[ps], 'M3', sep = '_')]] = patient_mids$SERc3 / patient_mids$PGc3
cortex_plasma[[paste(patient_site_names[ps], 'M1', sep = '_')]] = patient_mids$SERc1 / patient_mids$SERp1


}

p = ggplot(glioma_denovo, aes_string('Patient1E', 'Patient1E_M3')) +
  geom_point() +
  stat_poly_line(formula = y ~ x) +
  stat_poly_eq(use_label(c("eq", "adj.R2", "p.value")), formula = y ~ x, 
               label.x = "right") +
  ggtitle(paste0("spearman's correlation = ", 
                 round(cor(glioma_denovo$Patient1E, 
                           y = glioma_denovo$Patient1E_M3, method = 'pearson'), 
                       digits = 4))) +
  theme(text = element_text(size = 16))

p
ggsave(filename = paste('scatter', 'patient1E',
                        '.pdf', sep = '_'), 
       plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')

glioma_denovo_cor = calculateCorrelation(glioma_denovo)
pdf(file = 'corrplot_glioma_denovo.pdf', width = 10, height = 10)
p = correlationPlot(glioma_denovo_cor$M, glioma_denovo_cor$adj_pval)
print(p)
dev.off()

glioma_plasma_cor = calculateCorrelation(glioma_plasma)
pdf(file = 'corrplot_glioma_plasma.pdf', width = 10, height = 10)
p = correlationPlot(glioma_plasma_cor$M, glioma_plasma_cor$adj_pval)
print(p)
dev.off()

cortex_denovo_cor = calculateCorrelation(cortex_denovo)
pdf(file = 'corrplot_cortex_denovo.pdf', width = 10, height = 10)
p = correlationPlot(cortex_denovo_cor$M, cortex_denovo_cor$adj_pval)
print(p)
dev.off()

cortex_plasma_cor = calculateCorrelation(cortex_plasma)
pdf(file = 'corrplot_cortex_plasma.pdf', width = 10, height = 10)
p = correlationPlot(cortex_plasma_cor$M, cortex_plasma_cor$adj_pval)
print(p)
dev.off()

glioma_tme_cor = calculateCorrelation(glioma_tme)
pdf(file = 'corrplot_glioma_tme.pdf', width = 20, height = 20)
p = correlationPlot(glioma_tme_cor$M, glioma_tme_cor$adj_pval)
print(p)
dev.off()

glioma_denovo_ps = colMeans(glioma_denovo)
glioma_plasma_ps = colMeans(glioma_plasma)

cortex_denovo_ps = colMeans(cortex_denovo)
cortex_plasma_ps = colMeans(cortex_plasma)

glioma_ps = rbind(glioma_denovo_ps[16:30], 
                  glioma_plasma_ps[16:30],
                  cortex_denovo_ps[16:30],
                  cortex_plasma_ps[16:30])

glioma_ps_ML = rbind(glioma_denovo_ps[1:15],
                    glioma_plasma_ps[1:15],
                    cortex_denovo_ps[1:15],
                    cortex_plasma_ps[1:15])

cortex_ps_ML = rbind(cortex_denovo_ps,
                     cortex_plasma_ps)

rownames(cortex_ps_ML) = c('denovo_c', 'plasma_c')

rownames(glioma_ps) = c('SERg3/PGg3', 
                        'SERg1/SERp1', 
                        'SERc3/PGc3', 
                        'SERc1/SERp1')

rownames(glioma_ps_ML) = c('denovo_g', 
                           'plasma_g', 
                           'denovo_c', 
                           'plasma_c')

scaled_glioma_ps = scale(t(glioma_ps))
scaled_glioma_ps_ML = scale(t(glioma_ps_ML))
scaled_cortex_ps_ML = scale(t(cortex_ps_ML))
  
pdf(file = 'pheatmap_patient_MID_ratio.pdf', width = 10, height = 3)
pheatmap(t(scaled_glioma_ps), cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = COL1('OrRd'))
dev.off()

pdf(file = 'pheatmap_patient_ML_ratio_OrRd.pdf', width = 10, height = 3)
pheatmap(t(scaled_glioma_ps_ML), cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = COL1('OrRd'))#YlGn
dev.off()

cortex_denovo = read.delim(file = paste(ML_pred_dir, 'ML_denovo_cortex_1615_NandE.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
cortex_plasma = read.delim(file = paste(ML_pred_dir, 'ML_plasma_cortex_1615_NandE.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')

pdf(file = 'pheatmap_patient_ML_ratio_cortex_OrRd.pdf', width = 5, height = 3)
pheatmap(t(scaled_cortex_ps_ML), cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = COL1('OrRd'))
dev.off()

plasma_all = read.delim(file = paste(ML_pred_dir, 'ML_plasma_1327_2102_1615_NandE.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
scaled_denovo = scale(t(colMeans(denovo_all)))
pheatmap(t(scaled_denovo), cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = COL1('YlGn'))

#240523 - scatter plots - correlation of scaled flux & MID ratio
scaled_glioma_ps = data.frame(scaled_glioma_ps)
scaled_glioma_ps_ML = data.frame(scaled_glioma_ps_ML)

denovog_ser3 = data.frame('SERg3.PG3' = scaled_glioma_ps$SERg3.PGg3,
                          'denovo_g' = scaled_glioma_ps_ML$denovo_g,
                          'patient' = factor(seq(1, nrow(scaled_glioma_ps))))

plasmag_ser1 = data.frame('SERg1.SERp1' = scaled_glioma_ps$SERg1.SERp1,
                          'plasma_g' = scaled_glioma_ps_ML$plasma_g,
                          'patient' = factor(seq(1, nrow(scaled_glioma_ps))))
library("ggpubr")
p = ggscatter(denovog_ser3, x = "SERg3.PG3", y = "denovo_g", color = "patient",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "SERg3/PGg3", ylab = "Glucose-derived serine synthesis flux",
          label = "patient") +
  geom_smooth(method = "lm", color = "black") 
  
ggsave(filename = paste('scatter', 'denovo_c_SERc3_PGc3_cor',
                        '.pdf', sep = '_'), 
       plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')

p = ggscatter(plasmag_ser1, x = "SERg1.SERp1", y = "plasma_g", color = "patient",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "SERg1/SERp1", ylab = "Plasma serine uptake flux",
              label = "patient") +
  geom_smooth(method = "lm", color = "black") 

ggsave(filename = paste('scatter', 'plasma_c_SERc1_SERp1_cor',
                        '.pdf', sep = '_'), 
       plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')

# add cortex MIDs from nonenhancing and enhancing 
cortex_denovo = read.delim(file = paste(ML_pred_dir, 'ML_denovo_cortex_1615_NandE.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
cortex_plasma = read.delim(file = paste(ML_pred_dir, 'ML_plasma_cortex_1615_NandE.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')

for (ps in 1:(ncol(cortex_denovo)-1)){
  patient_mids = read.delim(file = paste(patient_folders[2*ps-1], '/mid_mc.txt', sep = ''), 
                            header = F, sep = '\t')
  
  colnames(patient_mids) = mid_names$V1
  E3 = patient_mids$SERc3 / patient_mids$PGc3
  E1 = patient_mids$SERc1 / patient_mids$SERp1
  
  patient_mids = read.delim(file = paste(patient_folders[2*ps], '/mid_mc.txt', sep = ''), 
                            header = F, sep = '\t')
  
  colnames(patient_mids) = mid_names$V1
  N3 = patient_mids$SERc3 / patient_mids$PGc3
  N1 = patient_mids$SERc1 / patient_mids$SERp1
  
  cortex_denovo[[paste(patient_site_names[2*ps-1], 'M3', sep = '_')]] = c(E3, N3)
  cortex_plasma[[paste(patient_site_names[2*ps-1], 'M1', sep = '_')]] = c(E1, N1)
  
}
ps = 8
patient_mids = read.delim(file = paste(patient_folders[2*ps-1], '/mid_mc.txt', sep = ''), 
                          header = F, sep = '\t')

colnames(patient_mids) = mid_names$V1
E3 = patient_mids$SERc3 / patient_mids$PGc3
E1 = patient_mids$SERc1 / patient_mids$SERp1

cortex_denovo[[paste(patient_site_names[2*ps-1], 'M3', sep = '_')]] = c(E3, E3)
cortex_plasma[[paste(patient_site_names[2*ps-1], 'M1', sep = '_')]] = c(E1, E1)

library(RcmdrMisc)
cortex_denovo_cor = calculateCorrelation(cortex_denovo)
pdf(file = 'corrplot_cortex_denovo_200.pdf', width = 10, height = 10)
p = correlationPlot(cortex_denovo_cor$M, cortex_denovo_cor$adj_pval)
print(p)
dev.off()

cortex_plasma_cor = calculateCorrelation(cortex_plasma)
pdf(file = 'corrplot_cortex_plasma_200.pdf', width = 10, height = 10)
p = correlationPlot(cortex_plasma_cor$M, cortex_plasma_cor$adj_pval)
print(p)
dev.off()


cortex_denovo_ps = colMeans(cortex_denovo)
cortex_plasma_ps = colMeans(cortex_plasma)

cortex_ps = rbind(cortex_denovo_ps[9:16],
                  cortex_plasma_ps[9:16])

cortex_ps_ML = rbind(cortex_denovo_ps[1:8],
                     cortex_plasma_ps[1:8])

rownames(cortex_ps) = c('SERg3/PGg3', 
                        'SERg1/SERp1')

rownames(cortex_ps_ML) = c('denovo_g', 
                           'plasma_g')

scaled_cortex_ps = scale(t(cortex_ps))
scaled_cortex_ps_ML = scale(t(cortex_ps_ML))

scaled_glioma_ps = scaled_cortex_ps
scaled_glioma_ps_ML = scaled_cortex_ps_ML
############################################################################## 
correlationPlot = function(M, adj_pval){
  p = corrplot(M, type = "full", method = 'ellipse', col = rev(COL2('PiYG')),
               tl.col = "black", tl.cex = 1.5, p.mat = adj_pval, 
               sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.1,
               insig = 'label_sig', diag = F) 
  return(p)
}

calculateCorrelation = function(df){
  M = cor(df, method = 'pearson')
  M = M[rowSums(is.na(M)) != ncol(M) - 1, colSums(is.na(M)) != nrow(M) - 1]
  
  adj_pval = rcorr.adjust(df)
  adj_pval = adj_pval$P
  adj_pval = ifelse(adj_pval == "<.0001", 0.0001, as.numeric(adj_pval))
  adj_pval = adj_pval[rowSums(is.na(adj_pval)) != ncol(adj_pval), 
                      colSums(is.na(adj_pval)) != nrow(adj_pval)]
  
  return(list('M' = M, 'adj_pval' = adj_pval))
}

##############################################################################
# scMFA-ML flux correlation ----
# scMFA folder
scMFA_dir = "C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/GBM/scRNA_Anjali/down_stream_analysis/darmanis2017/all_cells_batch_match_MFA"
patient_scMFA = list.files(path = scMFA_dir,
                            recursive = F)


neo_flux = c() # scMFA flux
for (ps in 1:length(patient_scMFA)){
  mal_mean = read.delim(file = paste(scMFA_dir, '/',
                                     patient_scMFA[ps], '/',
                                     patient_scMFA[ps], ' mal_mean.txt', 
                                     sep = ''), sep = ' ',
                        quote = '', row.names = 1, header = T)  
  neo_flux = c(neo_flux, list(mal_mean[c(1, 2, 3), ]))
}

neo_flux_df = data.frame(neo_flux)
colnames(neo_flux_df) = c('P1E', 'P1N', 'P2E', 'P4E', 'P4N', 
                          'P6E', 'P6N', 'P7E', 'P7N', 'P8N')
rownames(neo_flux_df) = c('denovo', 'plasma', 'tme')


neo_flux_normalized = sweep(neo_flux_df, 2, colSums(neo_flux_df),`/`)
neo_flux_scaled = scale(t(neo_flux_normalized))
neo_flux_scaled = data.frame(neo_flux_scaled)
neo_flux_normalized = data.frame(t(neo_flux_normalized))


# ML flux
glioma_denovo[, c('Patient2N', 'Patient3E',
                  'Patient3N', 'Patient5E',
                  'Patient5N')] = list(NULL)
glioma_plasma[, c('Patient2N', 'Patient3E',
                  'Patient3N', 'Patient5E',
                  'Patient5N')] = list(NULL)
glioma_tme[, c('Patient2N', 'Patient3E',
                  'Patient3N', 'Patient5E',
                  'Patient5N')] = list(NULL)

colnames(glioma_denovo) = c('P1E', 'P1N', 'P2E', 'P4E', 'P4N', 
                          'P6E', 'P6N', 'P7E', 'P7N', 'P8N')

colnames(glioma_plasma) = c('P1E', 'P1N', 'P2E', 'P4E', 'P4N', 
                            'P6E', 'P6N', 'P7E', 'P7N', 'P8N')

colnames(glioma_tme) = c('P1E', 'P1N', 'P2E', 'P4E', 'P4N', 
                            'P6E', 'P6N', 'P7E', 'P7N', 'P8N')

glioma_denovo_ps = colMeans(glioma_denovo)
glioma_plasma_ps = colMeans(glioma_plasma)
glioma_tme_ps = colMeans(glioma_tme)

glioma_ps = rbind(glioma_denovo_ps, glioma_plasma_ps, glioma_tme_ps)
scaled_glioma_ps = scale(t(glioma_ps))

scaled_glioma_ps = data.frame(scaled_glioma_ps)
colnames(scaled_glioma_ps) = c('denovo', 'plasma', 'tme')

glioma_ps = data.frame(t(glioma_ps))
colnames(glioma_ps) = c('denovo', 'plasma', 'tme')

neo_flux_df = data.frame(t(neo_flux_df))

# normalized
denovo_g = data.frame('denovo_ML' = glioma_ps$denovo,
                      'denovo_scMFA' = neo_flux_normalized$denovo,
                      'patient' = factor(rownames(glioma_ps)))

plasma_g = data.frame('plasma_ML' = glioma_ps$plasma,
                      'plasma_scMFA' = neo_flux_normalized$plasma,
                      'patient' = factor(rownames(glioma_ps)))

tme_g = data.frame('tme_ML' = glioma_ps$tme,
                   'tme_scMFA' = neo_flux_normalized$tme,
                   'patient' = factor(rownames(glioma_ps)))

# scaled
denovo_g = data.frame('denovo_ML' = scaled_glioma_ps$denovo,
                      'denovo_scMFA' = neo_flux_scaled$denovo,
                      'patient' = factor(rownames(scaled_glioma_ps)))

plasma_g = data.frame('plasma_ML' = scaled_glioma_ps$plasma,
                      'plasma_scMFA' = neo_flux_scaled$plasma,
                      'patient' = factor(rownames(scaled_glioma_ps)))

tme_g = data.frame('tme_ML' = scaled_glioma_ps$tme,
                   'tme_scMFA' = neo_flux_scaled$tme,
                   'patient' = factor(rownames(scaled_glioma_ps)))

library("ggpubr")
library(ggplot2)

p = ggscatter(denovo_g, x = "denovo_ML", y = "denovo_scMFA", color = "patient",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "denovo_ML", ylab = "denovo_scMFA",
              label = "patient") +
  geom_smooth(method = "lm", color = "black") 

ggsave(filename = 'scatter_denovo_ML_scMFA_glioma.pdf', 
       plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')

p = ggscatter(plasma_g, x = "plasma_ML", y = "plasma_scMFA", color = "patient",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "plasma_ML", ylab = "plasma_scMFA",
              label = "patient") +
  geom_smooth(method = "lm", color = "black") 

ggsave(filename = 'scatter_plasma_ML_scMFA_glioma.pdf', 
       plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')

p = ggscatter(tme_g, x = "tme_ML", y = "tme_scMFA", color = "patient",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "tme_ML", ylab = "tme_scMFA",
              label = "patient") +
  geom_smooth(method = "lm", color = "black") 

ggsave(filename = 'scatter_tme_ML_scMFA_glioma.pdf', 
       plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')

###############################################################################
# cBioPortal PanCancer ICGC/TCGA IDH-mutant vs IDH-wt glioma
shmt1_exp = read.delim('SHMT1_IDH1_mut_vs_wt.txt', header = T, sep = '\t')
shmt2_exp = read.delim('SHMT2_IDH1_mut_vs_wt.txt', header = T, sep = '\t')

colnames(shmt1_exp) = c('sample', 'IDH1_mut', 'exp', 'mutation')
colnames(shmt2_exp) = c('sample', 'IDH1_mut', 'exp', 'mutation')

shmt1_exp[['gene']] = 'SHMT1'
shmt2_exp[['gene']] = 'SHMT2'

shmt_exp = bind_rows(shmt1_exp, shmt2_exp)

stat_shmt1 = compare_means(data = shmt1_exp, 
                           formula = exp ~ IDH1_mut, 
                           method = 't.test',
                           p.adjust.method = 'bonf')

stat_shmt2 = compare_means(data = shmt2_exp, 
                           formula = exp ~ IDH1_mut, 
                           method = 't.test',
                           p.adjust.method = 'bonf')
stat_shmt1
stat_shmt2

p = ggplot(shmt1_exp, aes(x = gene, y = exp, fill = IDH1_mut)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.6)+
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5,
               position = position_dodge(0.8),
               binwidth = 1) +
  stat_summary(fun=mean, geom="point", size=1, shape = 4, color="black") +
  theme_classic()

ggsave(filename = 'TCGA_shmt1_IDH1.pdf', 
       plot = p, device = 'pdf', width = 4, height = 5, units = 'in')

###############################################################################
# PHGDH, PSAT1, PSPH gene expression and ratio of SERM+3/PGM+3

patient_folders = list.dirs(path = "C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/data_simulation/patient_MFA/patient_serine_MFA/Output files/serine_MFA_second_try_scflux",
                            recursive = F)
mid_names = read.delim(file = "C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/data_simulation/patient_MFA/patient_serine_MFA/mid_name.txt", 
                       header = F, sep = '\t')

patient_sites = c('P1E', 'P1N', 'P2E', 'P2N', 'P3E', 'P3N', 'P4E', 'P4N',
                  'P5E', 'P5N', 'P6E', 'P6N', 'P7E', 'P7N', 'P8N')

ratio_M3 = data.frame(matrix(nrow = 100, ncol = 0))

for (ps in 1:length(patient_folders)){
  patient_mids = read.delim(file = paste(patient_folders[ps], 
                                         '/mid_mc.txt', sep = ''), 
                            header = F, sep = '\t')
  
  
  colnames(patient_mids) = mid_names$V1
  
  ratio_M3[[patient_sites[ps]]] = patient_mids$SER3 / patient_mids$PG3
  
}

mean_M3 = colMeans(ratio_M3)
scaled_mean_M3 = scale(mean_M3)

available_scrna = scaled_mean_M3[c('P1E', 'P1N', 'P2E', 'P4E', 'P4N', 
                                   'P6E', 'P6N', 'P7E', 'P7N', 'P8N'), ]
# from darmanis.R
p = DotPlot(gbm, features = c('PSPH'), idents = c('Neoplastic'), 
            group.by = 'patient_site')

scdata = data.frame(patient_scid = as.character(p$data$id), gene_exp = p$data$avg.exp.scaled)
scdata[scdata$patient_scid == '6473SS_NonEnhancing', 'patient_scid'] = 'P8N'
scdata[scdata$patient_scid == '6369SS_NonEnhancing', 'patient_scid'] = 'P7N'
scdata[scdata$patient_scid == '6369SS_Enhancing', 'patient_scid'] = 'P7E'
scdata[scdata$patient_scid == '6025SS_NonEnhancing', 'patient_scid'] = 'P6N'
scdata[scdata$patient_scid == '6025SS_Enhancing', 'patient_scid'] = 'P6E'
scdata[scdata$patient_scid == '5675SS_NonEnhancing', 'patient_scid'] = 'P4N'
scdata[scdata$patient_scid == '5675SS_Enhancing', 'patient_scid'] = 'P4E'
scdata[scdata$patient_scid == '5363SR_Enhancing', 'patient_scid'] = 'P2E'
scdata[scdata$patient_scid == '4554DH_NonEnhancing', 'patient_scid'] = 'P1N'
scdata[scdata$patient_scid == '4554DH_Enhancing', 'patient_scid'] = 'P1E'
rownames(scdata) = scdata$patient_scid
sorted_scdata = scdata[names(available_scrna),]

gene_M3 = data.frame(gene_exp = sorted_scdata$gene_exp, 
                     ratio_M3 = available_scrna,
                     patient = names(available_scrna),
                     row.names = names(available_scrna))

library("ggpubr")
p = ggscatter(gene_M3, x = "gene_exp", y = "ratio_M3", color = "patient",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "PSPH", ylab = "SERg3/PGg3",
              label = "patient") +
  geom_smooth(method = "lm", color = "black") 

ggsave(filename = 'scatter_PSPH_serM3_PG3.pdf',
       plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')
