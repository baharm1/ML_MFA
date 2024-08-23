# Prepare input data for serine 13C-scMFA

# Import Libraries ------------------------------------------------------------
library(Seurat)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(Rmagic)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(cowplot)
library(tidyr)
library(patchwork)
library(reshape2)
library(Polychrome)  
library(scales)

##############################################################################
# heatmap of scfluxes to compare with ML fluxes
patient_scMFA = list.files('all_cells_batch_match_MFA')
neo_flux = c()
for (ps in 1:length(patient_scMFA)){
  mal_mean = read.delim(file = paste('all_cells_batch_match_MFA/', 
                                     patient_scMFA[ps], '/',
                                     patient_scMFA[ps], ' mal_mean.txt', 
                                     sep = ''), sep = ' ',
                        quote = '', row.names = 1, header = T)  
  neo_flux = c(neo_flux, list(mal_mean[c(1, 2, 3), ]))
}

neo_flux_df = data.frame(neo_flux)
colnames(neo_flux_df) = paste(rep('P', 10), seq(1, 10), sep = '')
rownames(neo_flux_df) = c('PG_SER', 'SERp_SER', 'SER_tran')


neo_flux_scaled = sweep(neo_flux_df, 2, colSums(neo_flux_df),`/`)

pdf(file = 'pheatmap_scMFA_scaled_compare_ML_Blues.pdf', width = 5, height = 3)
pheatmap(neo_flux_scaled, cluster_rows = F, cluster_cols = F, 
         show_rownames = T,
         color = COL1('Blues'))
dev.off()

###############################################################################
p = VlnPlot(gbm, features = c('ULK1', 'ATG13', 'ATG101', 'RB1CC1'), ncol = 2, pt.size = 0, group.by = 'cell_anno_2')
ggsave(filename = 'ULK1_complex_darmanis.pdf', plot = p, device = 'pdf', width = 8, height = 6)
p = VlnPlot(seuobj, features = c('ULK1', 'ATG13', 'C12orf44', 'RB1CC1'), ncol = 2, pt.size = 0)
###############################################################################
# save input scRNA patient data for purine metabolism
list_purine_genes = c('PPAT', 'GART', 'PFAS', 'PAICS', 'ADSL', 'ATIC', 
                      'PRPS1L1', 'PRPS1', 'PRPS2', # de novo IMP
                      'HPRT1', # HPX_IMP
                      'NT5C2', 'NT5C', 'NT5E', 'NT5DC4', 'NT5M', 
                      'NT5C1A', 'NT5C1B', # IMP_out
                      'ADSS', 'ADSSL1', # IMP_AMP 
                      'AMPD1', 'AMPD2', 'AMPD3', # AMP_IMP
                      'ADK', # ADE_AMP
                      'ADA', 'ADA2', 'LACC1', # ADE_INO
                      'PNP', # nucleoside to nucleobase
                      'IMPDH1', 'IMPDH2', 'GMPS', # IMP_GMP
                      'GUK1', # GMP_GDP
                      'RRM2B', 'RRM1', 'RRM2', 'NME6', 'AK9', 'NME7', 'NME1', 
                      'NME2', 'NME3', 'NME4', # GDP_out
                      'AK1', 'AK2', 'AK3', 'AK4', 'AK5', 'AK6', 'AK7', 'AK8', # AMP_out
                      'SLC28A1', 'SLC28A2', 'SLC28A3', # CNT
                      'SLC29A1', 'SLC29A2', 'SLC29A3', # ENT
                      'SLC43A3') # ENBT1
                      
pdf(file = 'dotplot_scwna_purine_input_scFBA_MFA_2406.pdf', height = 8.5, width = 5.5)
Clustered_DotPlot(gbm, features = list_purine_genes, group.by = 'cell_anno_2',
                  colors_use_exp = 
                    colorRampPalette(rev(brewer.pal(11, name = 'RdYlBu')))(12))
dev.off()

gbm = SetIdent(gbm, value = 'patientSite')
for (patient_site in unique(gbm$patientSite)){ 
  gbm_ps = subset(gbm, idents = patient_site)
  p = DotPlot_scCustom(gbm_ps, features = list_purine_genes, group.by = 'cell_anno_2',
                       flip_axes = T, x_lab_rotate = T) +
    theme(text = element_text(size = 14)) +
    scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                           limits = c(-2.5, 2.5))
  ggsave(filename = paste('dotplot_gbm_purine_list', patient_site, '.pdf', sep = '_'),
         plot = p, device = 'pdf', width = 6, height = 9)
}

gbm@meta.data[['patientSite']] = paste(gbm$patient, gbm$site, sep = '_')
table(gbm$patientSite)
mal = subset(gbm, idents = c('AC-like', 'NPC-like', 'MES-like', 'OPC-like', 'CellCycle'))
mal = SetIdent(mal, value = 'patientSite')

myeloid = subset(gbm, idents = c('Myeloid'))
myeloid = SetIdent(myeloid, value = 'patientSite')

rm(gbm)
gc()
for (patient_site in unique(mal$patientSite)){ 
  
  mal_ps = subset(mal, idents = patient_site)
  mal_ge = GetAssayData(mal_ps, slot = 'data')
  mal_ge = mal_ge[list_purine_genes, ]
  dim(mal_ge) 
  dimnames(mal_ge)[[1]] = paste(dimnames(mal_ge)[[1]], 'g', sep = '_')
  mal_ge = as.data.frame(mal_ge)
  
  my_ps = subset(myeloid, idents = patient_site)
  my_ge = GetAssayData(my_ps, slot = 'data')
  my_ge = my_ge[list_purine_genes, ]
  dim(my_ge) 
  dimnames(my_ge)[[1]] = paste(dimnames(my_ge)[[1]], 'm', sep = '_')
  my_ge = as.data.frame(my_ge)
  
  combined = bind_rows(my_ge, mal_ge)
  combined[is.na(combined)] = 0
  
  myeloid_ids = colnames(my_ps)
  neoplastic_ids = colnames(mal_ps)
  ids = c(myeloid_ids, neoplastic_ids)
  
  #cellnames_types = data.frame(ids, row.names = ids)
  #cellnames_types[['cell_type']] = NA
  
  #cellnames_types[cellnames_types$ids %in% myeloid_ids, 'cell_type'] = 'Myeloid'
  #cellnames_types[cellnames_types$ids %in% neoplastic_ids, 'cell_type'] = 'Neoplastic'
  #write.table(cellnames_types, file = paste(patient_site, 'purine_cellnames_types.txt', 
  #                                          sep = '_'), sep = '\t', quote = F)
  
  write.csv(combined, file = paste(patient_site, 'purine_scRNA.csv', sep = '_'), 
            quote = F)
  
}

###############################################################################
###############################################################################
# visualize output of scFBA-MFA for purine fluxes of neoplastic cells
flux_file = 'P1E_purine_scRNA_fg_20240628-150659'
flux = read.csv(file = paste('C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/Anjali/scfea/integrated_FBA_MFA/output_purine/', flux_file, '.csv', sep = ''), 
                row.names = 1)

flux = as.matrix(flux)
pheatmap(flux, cluster_rows = F, cluster_cols = F, show_rownames = F)

#median_neo_fluxes = apply(flux, 2, median)
mean_neo_fluxes = apply(flux, 2, mean)
sd_neo_fluxes = apply(flux, 2, sd)

#write.table(median_neo_fluxes, paste(flux_file, 'mal_purine_median.txt'))
write.table(mean_neo_fluxes, paste(flux_file, 'mal_purine_mean.txt'))
write.table(sd_neo_fluxes, paste(flux_file, 'mal_purine_sd.txt'))

cellnames_types = 'P1E_purine_cellnames_types'
cellnames_types = read.table(file = paste('C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/Anjali/scfea/integrated_FBA_MFA/patient_input_purine/patient_celltypes/', cellnames_types, '.txt', sep = ''), 
                             row.names = 1)
myeloid_ids = cellnames_types[cellnames_types$cell_type == 'Myeloid', 'ids']
neoplastic_ids = cellnames_types[cellnames_types$cell_type == 'Neoplastic', 'ids']

flux[['cell_type']] = NA
flux[myeloid_ids, 'cell_type'] = 'Myeloid'
flux[neoplastic_ids, 'cell_type'] = 'Neoplastic'

combine_myeloid_neo_f = function(f, new_str_f, old_str_f_m, old_str_f_n){
  f[[new_str_f]] = NA
  f[myeloid_ids, new_str_f] = f[myeloid_ids, old_str_f_m]
  f[neoplastic_ids, new_str_f] = f[neoplastic_ids, old_str_f_n]
  return(f)
}

flux = combine_myeloid_neo_f(flux, 'denovo_IMP', 'denovo_IMP_m', 'denovo_IMP_g')
flux = combine_myeloid_neo_f(flux, 'HPX_IMP', 'HPX_IMP_m', 'HPX_IMP_g')
flux = combine_myeloid_neo_f(flux, 'IMP_INO', 'IMP_INO_m', 'IMP_INO_g')
flux = combine_myeloid_neo_f(flux, 'IMP_AMP', 'IMP_AMP_m', 'IMP_AMP_g')
flux = combine_myeloid_neo_f(flux, 'AMP_IMP', 'AMP_IMP_m', 'AMP_IMP_g')
flux = combine_myeloid_neo_f(flux, 'IMP_GMP', 'IMP_GMP_m', 'IMP_GMP_g')
flux = combine_myeloid_neo_f(flux, 'ADE_AMP', 'ADE_AMP_m', 'ADE_AMP_g')
flux = combine_myeloid_neo_f(flux, 'AMP_ADE', 'AMP_ADE_m', 'AMP_ADE_g')
flux = combine_myeloid_neo_f(flux, 'AMP_out', 'AMP_out_m', 'AMP_out_g')
flux = combine_myeloid_neo_f(flux, 'ADE_in', 'ADE_in_m', 'ADE_in_g')
flux = combine_myeloid_neo_f(flux, 'INO_in', 'INO_in_m', 'INO_in_g')
flux = combine_myeloid_neo_f(flux, 'ADE_INO', 'ADE_INO_m', 'ADE_INO_g')
flux = combine_myeloid_neo_f(flux, 'INO_HPX', 'INO_HPX_m', 'INO_HPX_g')
flux = combine_myeloid_neo_f(flux, 'HPX_INO', 'HPX_INO_m', 'HPX_INO_g')
flux = combine_myeloid_neo_f(flux, 'HPX_in', 'HPX_in_m', 'HPX_in_g')
flux = combine_myeloid_neo_f(flux, 'GMP_GDP', 'GMP_GDP_m', 'GMP_GDP_g')
flux = combine_myeloid_neo_f(flux, 'GMP_GUO', 'GMP_GUO_m', 'GMP_GUO_g')
flux = combine_myeloid_neo_f(flux, 'GUA_GMP', 'GUA_GMP_m', 'GUA_GMP_g')
flux = combine_myeloid_neo_f(flux, 'GUO_in', 'GUO_in_m', 'GUO_in_g')
flux = combine_myeloid_neo_f(flux, 'GUO_GUA', 'GUO_GUA_m', 'GUO_GUA_g')
flux = combine_myeloid_neo_f(flux, 'GUA_GUO', 'GUA_GUO_m', 'GUA_GUO_g')
flux = combine_myeloid_neo_f(flux, 'GUA_in', 'GUA_in_m', 'GUA_in_g')
flux = combine_myeloid_neo_f(flux, 'GDP_out', 'GDP_out_m', 'GDP_out_g')


flux$cell_type = factor(flux$cell_type, levels = c('Myeloid', 'Neoplastic'))

plot_vln_myloid_neo = function(plt_str){
  p = ggplot(flux, aes_string(x = 'cell_type', y = plt_str, fill = 'cell_type')) +
    geom_violin(trim = F, scale = 'width') +
    geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
    stat_summary(fun=mean, geom="point", size=1, color="#FF0000") +
    stat_compare_means(label.y = c(0.01), size = 2) +
    
    #stat_compare_means(method = 't.test', size = 2) +
    theme_classic(base_size = 14) +
    theme(legend.position="none") +
    scale_fill_manual(values=c("#e7f582", "#f28500"))
  
  ggsave(filename = paste(flux_file, 'vlnplot', plt_str, '.pdf', 
                          sep = '_'), 
         plot = p, device = 'pdf', width = 2.5, height = 3)
  
}
plot_vln_myloid_neo('denovo_IMP')
plot_vln_myloid_neo('HPX_IMP')
plot_vln_myloid_neo('IMP_INO')
plot_vln_myloid_neo('IMP_AMP')
plot_vln_myloid_neo('AMP_IMP')
plot_vln_myloid_neo('IMP_GMP')
plot_vln_myloid_neo('ADE_AMP')
plot_vln_myloid_neo('AMP_out')
plot_vln_myloid_neo('ADE_in')
plot_vln_myloid_neo('INO_in')
plot_vln_myloid_neo('ADE_INO')
plot_vln_myloid_neo('INO_HPX')
plot_vln_myloid_neo('HPX_INO')
plot_vln_myloid_neo('HPX_in')
plot_vln_myloid_neo('GMP_GDP')
plot_vln_myloid_neo('GMP_GUO')
plot_vln_myloid_neo('GUA_GMP')
plot_vln_myloid_neo('GUO_in')
plot_vln_myloid_neo('GUO_GUA')
plot_vln_myloid_neo('GUA_GUO')
plot_vln_myloid_neo('GUA_in')
plot_vln_myloid_neo('GDP_out')


sub_in_flux = flux[flux$cell_type == 'Myeloid', c('denovo_IMP',
                                                  'HPX_IMP',
                                                  'IMP_INO',
                                                  'IMP_AMP',
                                                  'AMP_IMP',
                                                  'IMP_GMP',
                                                  'ADE_AMP',
                                                  'AMP_out',
                                                  'ADE_in',
                                                  'INO_in',
                                                  'ADE_INO',
                                                  'INO_HPX',
                                                  'HPX_in',
                                                  'GMP_GDP',
                                                  'GMP_GUO',
                                                  'GUA_GMP',
                                                  'GUO_in',
                                                  'GUO_GUA',
                                                  'GUA_in',
                                                  'GDP_out')]

sub_in_long = tidyr::pivot_longer(sub_in_flux, cols = 1:20)
sub_in_long[['col']] = factor(sub_in_long$name, 
                              levels = c('denovo_IMP',
                                         'HPX_IMP',
                                         'IMP_INO',
                                         'IMP_AMP',
                                         'AMP_IMP',
                                         'IMP_GMP',
                                         'ADE_AMP',
                                         'AMP_out',
                                         'ADE_in',
                                         'INO_in',
                                         'ADE_INO',
                                         'INO_HPX',
                                         'HPX_in',
                                         'GMP_GDP',
                                         'GMP_GUO',
                                         'GUA_GMP',
                                         'GUO_in',
                                         'GUO_GUA',
                                         'GUA_in',
                                         'GDP_out'))

my_comparisons = list(c("denovo_IMP", "HPX_IMP"), 
                      c("IMP_AMP", "ADE_AMP"),
                      c("IMP_GMP", "GUA_GMP"))

p = ggplot(sub_in_long, aes_string(x = 'col', y = 'value', fill = 'col')) +
  geom_violin(trim = F, scale = 'width') +
  geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
  #ggtitle('PG_SER') +
  #stat_compare_means(comparisons = my_comparisons, 
                  #   label.y = c(0.1, 0.13, 0.16), size = 2) +
  #stat_summary(fun=mean, geom="point", size=1, color="#E69F00") +
  #stat_compare_means(method = 't.test', size = 2) +
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.2)) 
  
scale_fill_manual(values=c("#f36aa3", "#28A8E0", "#28A8E0", 
                             "#f36aa3", "#f36aa3", "#28A8E0", 
                             "#f36aa3", 
                             "#f36aa3", "#28A8E0", "#f36aa3", 
                             "#f36aa3", "#f36aa3"))


ggsave(filename = paste(flux_file, 'vlnplot', 'myeloid_purine_fluxes', '.pdf', 
                        sep = '_'), 
       plot = p, device = 'pdf', width = 5, height = 3.5)

# half violin plots
sub_in_flux = flux[, c('cell_type', 
                       'denovo_IMP',
                       'HPX_IMP',
                       'IMP_INO',
                       'IMP_AMP',
                       'AMP_IMP',
                       'IMP_GMP',
                       'ADE_AMP',
                       'AMP_ADE',
                       'AMP_out',
                       'ADE_in',
                       'INO_in',
                       'ADE_INO',
                       'INO_HPX',
                       'HPX_INO',
                       'HPX_in',
                       'GMP_GDP',
                       'GMP_GUO',
                       'GUA_GMP',
                       'GUO_in',
                       'GUO_GUA',
                       'GUA_GUO',
                       'GUA_in',
                       'GDP_out')]
sub_in_long = tidyr::pivot_longer(sub_in_flux, cols = 2:24)
sub_in_long$name = factor(sub_in_long$name, 
                              levels = c('denovo_IMP',
                                         'HPX_IMP',
                                         'IMP_INO',
                                         'IMP_AMP',
                                         'AMP_IMP',
                                         'IMP_GMP',
                                         'ADE_AMP',
                                         'AMP_ADE',
                                         'AMP_out',
                                         'ADE_in',
                                         'INO_in',
                                         'ADE_INO',
                                         'INO_HPX',
                                         'HPX_INO',
                                         'HPX_in',
                                         'GMP_GDP',
                                         'GMP_GUO',
                                         'GUA_GMP',
                                         'GUO_in',
                                         'GUO_GUA',
                                         'GUA_GUO',
                                         'GUA_in',
                                         'GDP_out'))

p = ggplot(sub_in_long, aes(x = name, y = value, fill = cell_type)) +
  introdataviz::geom_split_violin(alpha = 0.7, scale = 'width') +
  geom_boxplot(outlier.shape = NA, width = 0.3,
               show.legend = FALSE) +
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = 'red') +
  scale_x_discrete(name = "Purine Fluxes") +
  scale_fill_manual(values = c("#e7f582", "#f28500"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

ggsave(filename = paste(flux_file, 
                        'vlnplot_myeloid_neoplastic_purine_fluxes.pdf',
                        sep = '_'), 
       plot = p, device = 'pdf', height = 4, width = 8, units = 'in')



sub_in_flux = flux[flux$cell_type == 'Neoplastic', c('denovo_IMP',
                                                     'HPX_IMP',
                                                     'IMP_out',
                                                     'IMP_AMP',
                                                     'AMP_IMP',
                                                     'ADE_AMP',
                                                     'IMP_GMP',
                                                     'GMP_GDP',
                                                     'GMP_out',
                                                     'GUA_GMP',
                                                     'GDP_out',
                                                     'AMP_out')]

sub_in_long = pivot_longer(sub_in_flux, cols = 1:12)
sub_in_long[['col']] = factor(sub_in_long$name, 
                              levels = c('denovo_IMP', 'HPX_IMP', 'AMP_IMP',
                                         'IMP_out', 'IMP_AMP','ADE_AMP',
                                         'AMP_out',
                                         'IMP_GMP', 'GUA_GMP', 'GMP_GDP', 
                                         'GMP_out', 'GDP_out'))

p = ggplot(sub_in_long, aes_string(x = 'col', y = 'value', fill = 'col')) +
  geom_violin(trim = F, scale = 'width') +
  geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
  #ggtitle('PG_SER') +
  stat_compare_means(comparisons = my_comparisons, 
                     label.y = c(0.1, 0.13, 0.16), size = 2) +
  stat_summary(fun=mean, geom="point", size=1, color="#E69F00") +
  #stat_compare_means(method = 't.test', size = 2) +
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.2)) +
  scale_fill_manual(values=c("#f36aa3", "#28A8E0", "#28A8E0", 
                             "#f36aa3", "#f36aa3", "#28A8E0", 
                             "#f36aa3", 
                             "#f36aa3", "#28A8E0", "#f36aa3", 
                             "#f36aa3", "#f36aa3"))


ggsave(filename = paste(flux_file, 'vlnplot', 'neoplastic_purine_fluxes', '.pdf', 
                        sep = '_'), 
       plot = p, device = 'pdf', width = 5, height = 3.5)

###############################################################################
# expression of nucleoside transporters

list_purine_genes = c('PPAT', 'GART', 'PFAS', 'PAICS', 'ADSL', 'ATIC', 
                      'PRPS1L1', 'PRPS1', 'PRPS2', # de novo IMP
                      'HPRT1', # HPX_IMP
                      'LACC1', 'PNP', # INO_HPX and GUO_GUA
                      'NT5C2', 'NT5C', 'NT5E', 'NT5DC4', 'NT5M', 
                      'NT5C1A', 'NT5C1B', # IMP_out
                      'ADSS', 'ADSSL1', # IMP_AMP 
                      'AMPD1', 'AMPD2', 'AMPD3', # AMP_IMP
                      'ADK', # ADE_AMP
                      'IMPDH1', 'IMPDH2', 'GMPS', # IMP_GMP
                      'GUK1', # GMP_GDP
                      'RRM2B', 'RRM1', 'RRM2', 'NME6', 'AK9', 'NME7', 'NME1', 
                      'NME2', 'NME3', 'NME4', # GDP_out
                      'AK1', 'AK2', 'AK3', 'AK4', 'AK5', 'AK6', 'AK7', 'AK8', # AMP_out
                      'SLC28A1', 'SLC28A2', 'SLC28A3',
                      'SLC29A1', 'SLC29A2', 'SLC29A3', 'SLC29A4')

gbm$cell_anno_2 = factor(gbm$cell_anno_2, 
                         levels = c('AC-like', 'NPC-like', 
                                    'OPC-like', 'MES-like',  'CellCycle',  
                                    'Oligo', 'Endothelial', 'Pericyte', 
                                    'Myeloid', 'Lymphoid'))

pdf(file = 'dotplot_scwna_purine_model_new.pdf', height = 11, width = 5.5)
DotPlot_scCustom(gbm, features = list_purine_genes, group.by = 'cell_anno_2',
                 flip_axes = T, x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))
dev.off()

pdf(file = 'dotplot_darmanis_purine_model_new.pdf', height = 11, width = 5.5)
DotPlot_scCustom(seuobj, features = list_purine_genes, group.by = 'major_cell_type',
                 flip_axes = T, x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))
dev.off()

p = DotPlot_scCustom(gbm, features = c("SLC28A1",
                                       "SLC28A2", 
                                       "SLC28A3",
                                       "SLC29A1", 
                                       "SLC29A2",
                                       "SLC29A3"), 
                     group.by = 'cell_anno_2',
                     flip_axes = T, x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))
ggsave(filename = 'dotplot_scwna_nucleoside_transporter.pdf', plot = p, 
       device = 'pdf', height = 4, width = 6, units = 'in')

p = DotPlot_scCustom(seuobj, features = c("SLC28A1",
                                          "SLC28A2", 
                                          "SLC28A3",
                                          "SLC29A1", 
                                          "SLC29A2",
                                          "SLC29A3"), 
                     group.by = 'major_cell_type',
                     flip_axes = T, x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))
ggsave(filename = 'dotplot_darmanis_nucleoside_transporter.pdf', plot = p, 
       device = 'pdf', height = 4, width = 6, units = 'in')


magic_obj = magic(t(GetAssayData(gbm, slot = 'data')), genes = c("SLC29A1", 
                                                                 "SLC29A2",
                                                                 "SLC29A3",
                                                                 "SLC29A4",
                                                                 "SLC28A1",
                                                                 "SLC28A2", 
                                                                 "SLC28A3"))
magic_res = magic_obj$result
###############################################################################
# gene expression of lactate and glutamine production and consumption for grant purposes

lactate_prod = c('LDHA', 'LDHAL6A', 'LDHAL6B', 'LDHB', 'LDHC', 'LDHD')

lactate_tran = c('SLC16A1', 'SLC16A3', 'SLC16A4', 'SLC5A8', 'SLC5A12')

gln_prod = c('GLUL' # Glu == Gln
             )
gln_cons = c('PPAT', #'GART', 'PFAS', 'PAICS', 'ADSL', 'ATIC', 
             #'PRPS1L1', 'PRPS1', 'PRPS2', # PRPP + Gln == IMP
             'CPS1', 'CAD', #'DHODH', 'UMPS', # GLn == UMP
             'GLS', 'GLS2', # Gln == Glu
             'GFPT1', # F6P == Glucosamine 6P
             'ASNS'
             )

gln_tran = c('SLC1A5', 'SLC3A2', 'SLC7A5', 'SLC6A14', 'SLC38A1', 'SLC38A2')

bcaa_tran = c('SLC3A2', 'SLC7A5', 'SLC16A9', 'SLC25A44')
bcka_tran = c('SLC16A1', 'SLC16A3')
bcka_prod = c('BCAT1', 'BCAT2')
bcka_cons = c('BCKDHA', 'BCKDHB', 'DLD', 'DBT')



magic_obj = magic(t(GetAssayData(seuobj, slot = 'data')))
magic_res = magic_obj$result

lactate_prod_magic = rowSums(magic_res[, lactate_prod])
lactate_tran_magic = rowSums(magic_res[, lactate_tran])
gln_prod_magic = magic_res[, gln_prod]
gln_cons_magic = rowSums(magic_res[, gln_cons])
gln_tran_magic = rowSums(magic_res[, gln_tran])

gln_prod_cons = gln_prod_magic - gln_cons_magic

bcaa_tran_magic = rowSums(magic_res[, bcaa_tran])
bcka_tran_magic = rowSums(magic_res[, bcka_tran])
bcka_prod_magic = rowSums(magic_res[, bcka_prod])
bcka_cons_magic = rowSums(magic_res[, bcka_cons])

bcka_prod_cons = bcaa_tran_magic + bcka_prod_magic - bcka_cons_magic
bcka_prod_cons2 = bcka_prod_magic - bcka_cons_magic

seuobj@meta.data[['lactate_prod_magic']] = lactate_prod_magic
seuobj@meta.data[['lactate_tran_magic']] = lactate_tran_magic
seuobj@meta.data[['gln_prod_magic']] = gln_prod_magic
seuobj@meta.data[['gln_cons_magic']] = gln_cons_magic
seuobj@meta.data[['gln_tran_magic']] = gln_tran_magic
seuobj@meta.data[['gln_prod_cons']] = gln_prod_cons

seuobj@meta.data[['bcaa_tran_magic']] = bcaa_tran_magic
seuobj@meta.data[['bcka_tran_magic']] = bcka_tran_magic
seuobj@meta.data[['bcka_prod_magic']] = bcka_prod_magic
seuobj@meta.data[['bcka_cons_magic']] = bcka_cons_magic
seuobj@meta.data[['bcka_prod_cons']] = bcka_prod_cons
seuobj@meta.data[['bcka_prod_cons2']] = bcka_prod_cons2

p = VlnPlot(seuobj, features = c('bcka_prod_cons2'), 
            group.by = 'major_cell_type',
            #cols = darmanis_colors,
) + 
  NoLegend() +
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white") #+
#geom_hline(yintercept = 0, linetype = "dashed", color = 'red')

ggsave(filename = 'vlnplot_major_cell_type_bcka_prod_cons2.pdf', plot = p, 
       device = 'pdf', height = 6, width = 6, units = 'in')

# 
p = VlnPlot(seuobj, features = c('serine_prod_cons'), 
            group.by = 'Location',
            #cols = darmanis_colors,
) + 
  NoLegend() +
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white") #+
#geom_hline(yintercept = 0, linetype = "dashed", color = 'red')

ggsave(filename = 'vlnplot_location_serine_prod_cons.pdf', plot = p, 
       device = 'pdf', height = 6, width = 3, units = 'in')


p = DotPlot_scCustom(seuobj, features = c(bcka_prod, bcka_cons, bcka_tran, bcaa_tran), group.by = 'major_cell_type',
                     flip_axes = T, x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_major_cell_type_bcka_prod_cons_tran.pdf', plot = p, 
       device = 'pdf', height = 4, width = 6, units = 'in')

###############################################################################
# serine transportes expression in our data

magic_obj = magic(t(GetAssayData(gbm, slot = 'data')), genes = serine_tran)
magic_res = magic_obj$result
magic_res[['cell_anno_2']] = gbm$cell_anno_2
magic_res[['cell_anno']] = as.character(magic_res$cell_anno_2)

magic_res[magic_res$cell_anno == 'AC-like', 'cell_anno'] = 'Neoplastic'
magic_res[magic_res$cell_anno == 'MES-like', 'cell_anno'] = 'Neoplastic'
magic_res[magic_res$cell_anno == 'OPC-like', 'cell_anno'] = 'Neoplastic'
magic_res[magic_res$cell_anno == 'NPC-like', 'cell_anno'] = 'Neoplastic'
magic_res[magic_res$cell_anno == 'CellCycle', 'cell_anno'] = 'Neoplastic'

my_colors = c('#E4E1E3',
              '#16FF32',
              '#3283FE',
              '#B00068',
              '#90AD1C',
              '#5A5156',
              '#FE00FA',
              '#FEAF16',
              '#1CFFCE',
              '#F6222E')

scales::show_col(my_colors)

my_colors_2 = c('#AAF400',
              '#5A5156',
              '#FE00FA',
              '#FEAF16',
              '#1CFFCE',
              '#F6222E')
magic_res$cell_anno = factor(magic_res$cell_anno, levels = c('Neoplastic', 'Myeloid',
                                                             'Oligo', 'Endothelial', 
                                                             'Pericyte', 'Lymphoid'))

my_comparisons = list(c("Neoplastic", "Endothelial"), c("Neoplastic", "Pericyte"))

  stat_compare_means(comparisons = my_comparisons, label.y = c(0.95, 0.75, 1), size = 2) +
  #stat_compare_means(method = 't.test', size = 2) +
  theme_classic(base_size = 14) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#338333", "#F28500", "#007BBB"))

plot_list = list()
i = 0
for (gene in serine_tran){
  i = i + 1
  
  p = ggplot(magic_res, aes_string(x = 'cell_anno', y = gene, fill = 'cell_anno')) +
    geom_violin(scale = 'width', trim = F) +
    geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
    #stat_compare_means(method = 't.test', size = 2) +
    theme_classic(base_size = 8) +
    theme(legend.position="none") +
    scale_fill_manual(values = my_colors_2)
  
  plot_list[[i]] = p
}
plot_gene = cowplot::plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
                      plot_list[[4]], plot_list[[5]], plot_list[[6]], 
                      plot_list[[7]], plot_list[[8]], plot_list[[9]], 
                      plot_list[[10]], plot_list[[11]], plot_list[[12]], 
                      plot_list[[13]], plot_list[[14]], plot_list[[15]], 
                      plot_list[[16]], plot_list[[17]],
                      nrow = 6, ncol = 3)
ggsave(filename = 'vlnplot_serine_tran_magic_major.pdf', 
       plot = plot_gene, device = 'pdf', width = 10, height = 15)

gbm@meta.data[['cell_anno']] = as.character(gbm$cell_anno_2)
gbm@meta.data[gbm$cell_anno == 'AC-like', 'cell_anno'] = 'Neoplastic'
gbm@meta.data[gbm$cell_anno == 'MES-like', 'cell_anno'] = 'Neoplastic'
gbm@meta.data[gbm$cell_anno == 'OPC-like', 'cell_anno'] = 'Neoplastic'
gbm@meta.data[gbm$cell_anno == 'NPC-like', 'cell_anno'] = 'Neoplastic'
gbm@meta.data[gbm$cell_anno == 'CellCycle', 'cell_anno'] = 'Neoplastic'

gbm$cell_anno = factor(gbm$cell_anno, levels = c('Neoplastic', 'Myeloid',
                                                 'Oligo', 'Endothelial', 
                                                 'Pericyte', 'Lymphoid'))

p = DotPlot_scCustom(gbm, features = serine_tran, group.by = 'cell_anno_2', 
                     idents = c('AC-like', 'NPC-like', 'OPC-like', 'MES-like', 'CellCycle'),
                     flip_axes = T, x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_serine_tran_neoplastic_scwna.pdf', plot = p, 
       device = 'pdf', height = 5, width = 5, units = 'in')
###############################################################################
# markers of seurat clusters
gbm = SetIdent(gbm, value = 'RNA_snn_res.0.1')
markers <- FindAllMarkers(gbm, only.pos = TRUE)
saveRDS(markers, file = 'markers_res_01.rds')

markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

pdf(file = 'heatmap_top20.pdf', width = 9, height = 60)
DoHeatmap(gbm, features = top20$gene) + 
  scale_fill_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu'))) 
                         #limits = c(-2.5, 2.5)) 
#  NoLegend()
dev.off()

markers_sig = markers[(markers$avg_log2FC > 1) & (markers$p_val_adj < 0.05), ]

markers_sig %>%
  group_by(cluster) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20_sig

pdf(file = 'heatmap_top20_2.pdf', width = 9, height = 20)
DoHeatmap(gbm, features = top20$gene) + 
  scale_fill_gradientn(colours = rev(brewer.pal(11, name = 'PiYG'))) 
#limits = c(-2.5, 2.5)) 
#  NoLegend()
dev.off()
###############################################################################
# pie plot cell type proportions in the whole dataset
prop = table(as.character(gbm$cell_anno_2))
no_all_cells = sum(prop)
prop = prop / no_all_cells * 100
prop = as.data.frame(prop)
rownames(prop) = prop$Var1
ordered_name = c('AC-like', 'NPC-like', 'OPC-like', 'MES-like', 'CellCycle',
                 'Myeloid', 'Oligo', 'Endothelial', 'Pericyte', 'Lymphoid')
prop = prop[ordered_name, ]

prop
pdf(file = 'pie_cell_anno_2.pdf', height = 3, width = 3)
pie(prop$Freq , labels = prop$Var1, border = "black", col = my_colors)
dev.off()

prop['AC-like', 'Freq'] + prop['NPC-like', 'Freq'] + prop['OPC-like', 'Freq'] + prop['MES-like', 'Freq'] + prop['CellCycle', 'Freq']
##############################################################################
# Sravya NAG project: myeloids secrete NAG for gliomas
seuobj = readRDS('darmanis.rds')
seuobj = SetIdent(seuobj, value = c('major_cell_type'))

sravya_genes = c('NAGS', 'SLC7A11', 'ACY1', 'ACSS2', 'CPS1',
                 'SLC1A1', 'SLC1A2', 'SLC1A3', 'SLC1A4',
                 'SLC1A6', 'SLC1A7', 'SLC17A5', 'SLC17A6', 'SLC17A7', 'SLC17A8',
                 'SLC25A12', 'SLC25A13', 'SLC25A22')
NAG_prod = c('NAGS')
NAG_cons = c('ACY1', 'ACSS2', 'CPS1')

NAG_tran = c('SLC1A1', 'SLC1A2', 'SLC1A3', 'SLC1A4',
             'SLC1A6', 'SLC1A7', 'SLC7A11',
             'SLC17A5', 'SLC17A6', 'SLC17A7', 'SLC17A8',
             'SLC25A12', 'SLC25A13', 'SLC25A22')

p = DotPlot_scCustom(subset_seuobj, features = c(NAG_prod, NAG_cons, NAG_tran), 
                     group.by = 'major_cell_type',
                     idents = c('Myeloid', 'Neoplastic', 'Oligo', 'OPC', 'Vascular'),
                     flip_axes = T, x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_sravya_tumor_all_genes_darmanis.pdf', plot = p, 
       device = 'pdf', height = 5, width = 5, units = 'in')


magic_obj = magic(t(GetAssayData(seuobj, slot = 'data')))
magic_res = magic_obj$result

NAG_prod_magic = rowSums(magic_res[, NAG_prod])
NAG_prod_magic = magic_res[, NAG_prod]
NAG_cons_magic = rowSums(magic_res[, NAG_cons])
NAG_tran_magic = rowSums(magic_res[, NAG_tran])

NAG_prod_cons = NAG_prod_magic - NAG_cons_magic

seuobj@meta.data[['NAG_prod_cons']] = NAG_prod_cons
seuobj@meta.data[['NAG_prod']] = NAG_prod_magic
seuobj@meta.data[['NAG_cons']] = NAG_cons_magic
seuobj@meta.data[['NAG_tran']] = NAG_tran_magic

seuobj = SetIdent(seuobj, value = 'Location')
subset_seuobj = subset(seuobj, idents = c('Tumor'))
subset_seuobj = SetIdent(subset_seuobj, value = 'major_cell_type')

p = VlnPlot(seuobj, features = c('NAG_cons'), 
            group.by = 'major_cell_type',
            #idents = c('Myeloid', 'Neoplastic', 'Oligo', 'OPC', 'Vascular'),
            #cols = darmanis_colors,
) + 
  NoLegend() +
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white") +
  
  stat_compare_means(comparisons = list(c('Neoplastic', 'Myeloid')), 
                     label.y = c(0), size = 2) +
  stat_summary(fun=mean, geom="point", size=1, color="#E69F00") 
#geom_hline(yintercept = 0, linetype = "dashed", color = 'red')

ggsave(filename = 'vlnplot_major_cell_type_NAG_cons.pdf', plot = p, 
       device = 'pdf', height = 6, width = 6, units = 'in')

# 
p = VlnPlot(seuobj, features = c('NAG_prod_cons'), 
            group.by = 'Location',
            #cols = darmanis_colors,
) + 
  NoLegend() +
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white") #+
#geom_hline(yintercept = 0, linetype = "dashed", color = 'red')

ggsave(filename = 'vlnplot_location_serine_prod_cons.pdf', plot = p, 
       device = 'pdf', height = 6, width = 3, units = 'in')


p = DotPlot_scCustom(seuobj, features = serine_cons, group.by = 'cell_type',
                     flip_axes = T, x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_cell_type_serine_cons.pdf', plot = p, 
       device = 'pdf', height = 5, width = 6, units = 'in')

###############################################################################
# MEBOCOST

commu_df = read.csv(file = 'C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/GBM/scRNA_Anjali/down_stream_analysis/revised_figs_2310/cellCommun/MEBOCOST/darmanis_tumor/commu_df.csv',
                    header = T, quote = '')

commu_df = commu_df[commu_df$permutation_test_fdr < 0.05, ]

metab_connectivity = as.data.frame(table(commu_df$Metabolite_Name))

metab_connectivity = metab_connectivity[order(metab_connectivity$Freq,
                                              decreasing = F),]
metab_connectivity$Var1 = factor(metab_connectivity$Var1, 
                                 levels = metab_connectivity$Var1)

p = ggplot(metab_connectivity, aes(x = Var1, y = Freq, fill = Var1))+
  geom_bar(width = 0.7, stat = "identity")+
  coord_flip() +
  theme_classic() +
  theme(legend.position="none") 
        #axis.text.x = element_text(angle = 45)) 
  
ggsave(filename = 'MEBOCOST_sig_metabolite_connectivity.pdf', plot = p, 
       device = 'pdf', height = 5, width = 4, units = 'in')

##############################################################################
# mice predictions
mice_denovo_g = read.delim(file = 'C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_glioma_denovo_4_20240408_1327_.txt',
                           header = T, quote = '', row.names = 1)
mice_denovo_g[['flux']] = 'denovo'
mice_plasma_g = read.delim(file = 'C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_glioma_plasma_4_20240408_2102_.txt',
                           header = T, quote = '', row.names = 1)
mice_plasma_g[['flux']] = 'plasma'
mice_denovo_c = read.delim(file = 'C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_cortex_denovo_4_20240508_1220_.txt',
                           header = T, quote = '', row.names = 1)
mice_denovo_c[['flux']] = 'denovo'

colnames(mice_denovo_g) = c('value', 'flux')
colnames(mice_plasma_g) = c('value', 'flux')
colnames(mice_denovo_c) = c('value', 'flux')

mice_g = rbind(mice_denovo_g, mice_plasma_g)

p = ggplot(mice_g, aes_string(x = 'flux', y = 'value', fill = 'flux')) +
  geom_violin(trim = F, scale = 'width') +
  geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
  #ggtitle('PG_SER') +
  #stat_compare_means(comparisons = my_comparisons, 
  #                   label.y = c(0.1, 0.13, 0.16), size = 2) +
  stat_summary(fun=mean, geom="point", size=1, color="#E69F00") +
  #stat_compare_means(method = 't.test', size = 2) +
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.2)) +
  scale_fill_manual(values=c("#f36aa3", "#28A8E0"))


ggsave(filename = 'mice_glioma_ML_pred.pdf',
       plot = p, device = 'pdf', width = 3, height = 3.5)

# mean and standard deviation of predicted fluxes of serine in mice
mice_denovo_g = read.delim(file = 'C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_glioma_denovo_4_20240408_1327_.txt',
                           header = T, quote = '', row.names = 1)
mice_plasma_g = read.delim(file = 'C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_glioma_plasma_4_20240408_2102_.txt',
                           header = T, quote = '', row.names = 1)
mice_denovo_c = read.delim(file = 'C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_cortex_denovo_4_20240508_1220_.txt',
                           header = T, quote = '', row.names = 1)

mice_g = cbind(mice_denovo_g, mice_plasma_g)
mice_g[['tme']] = 1 - mice_g$denovo - mice_g$plasma

mice_c = mice_denovo_c
mice_c[['plasma']] = 1 - mice_denovo_c$denovo

colMeans(mice_g)
colMeans(mice_c)

sd(mice_g$denovo)
sd(mice_g$plasma)
sd(mice_g$tme)

sd(mice_c$denovo)
sd(mice_c$plasma)

# MMF experiment
# mean and standard deviation of predicted fluxes of serine in mice
timestr = '20240510_1434'
mice_ctrl = read.delim(file = paste('C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_mmf__ctrl__4_',
                                    timestr, '_.txt', sep = ''),
                       header = T, quote = '', row.names = 1)
mice_1low = read.delim(file = paste('C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_mmf__Tx1_low__4_', 
                       timestr, '_.txt', sep = ''),
                       header = T, quote = '', row.names = 1)
mice_1high = read.delim(file = paste('C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_mmf__Tx1_high__4_', 
                                     timestr, '_.txt', sep = ''),
                           header = T, quote = '', row.names = 1)

mice_2low = read.delim(file = paste('C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_mmf__Tx2_low__4_', 
                                    timestr, '_.txt', sep = ''),
                       header = T, quote = '', row.names = 1)
mice_2high = read.delim(file = paste('C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_mmf__Tx2_high__4_', 
                                     timestr, '_.txt', sep = ''),
                        header = T, quote = '', row.names = 1)

mice_4low = read.delim(file = paste('C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_mmf__Tx4_low__4_', 
                                    timestr, '_.txt', sep = ''),
                       header = T, quote = '', row.names = 1)
mice_4high = read.delim(file = paste('C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_mice_mmf__Tx4_high__4_', 
                                     timestr, '_.txt', sep = ''),
                        header = T, quote = '', row.names = 1)

mice_g = cbind(mice_ctrl, mice_1low, mice_1high, 
               mice_2low, mice_2high, mice_4low, mice_4high)
colnames(mice_g) = c('ctrl', 'low_1', 'high_1', 'low_2', 'high_2', 'low_4', 'high_4')

colMeans(mice_g)
apply(mice_g, 2, sd)


# purine patient predictions

patient_dir = paste('C:/Users/bme/Dropbox (University of Michigan)/Lab_Baharan/ML_serine_model/pred_patient_purine_', 
                    timestr, sep = '')
patient_files = list.files(path = patient_dir)
patient_pred_all = data.frame(temp = rep(as.double(NA), 100))
for(ps in 1:length(patient_files)){
  patient_pred = read.delim(file = paste(patient_dir, '/', patient_files[ps], sep = ''), 
             header = T, quote = '', row.names = 1)
  colnames(patient_pred) = paste('P', ps, sep = '')
  patient_pred_all = cbind(patient_pred_all, patient_pred)
  
}

patient_pred_all$temp = NULL
colnames(patient_pred_all) = c('P1E', 'P1N', 
                               'P2E', 'P2N', 
                               'P3E', 'P3N', 
                               'P4E', 'P4N',
                               'P5E', 'P5N',
                               'P6E', 'P6N', 
                               'P7E', 'P7N', 
                               'P8N')

colMeans(patient_pred_all)
apply(patient_pred_all, 2, sd)
################################################################################
# Sravya - NAG ; Help Dan with grant figures
sravya_genes = c('NAGS', 'NAT8L', 'ASPA', 'ACY1', 'ACSS1', 'ACSS2')

p = DotPlot_scCustom(seuobj, features = sravya_genes, 
                     group.by = 'major_cell_type',
                     #idents = c('Myeloid', 'Neoplastic', 'Oligo', 'OPC', 'Vascular'),
                     #flip_axes = T, 
                     x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_darmanis_NAG_sravya.pdf', plot = p, 
       device = 'pdf', height = 4, width = 5, units = 'in')


df = data.frame(sample = gbm@assays$RNA@data@Dimnames[[2]], anno = gbm$cell_anno_2)
df[df$anno == 'AC-like', 'anno'] = 'Neoplastic'
df[df$anno == 'NPC-like', 'anno'] = 'Neoplastic'
df[df$anno == 'OPC-like', 'anno'] = 'Neoplastic'
df[df$anno == 'MES-like', 'anno'] = 'Neoplastic'
df[df$anno == 'CellCycle', 'anno'] = 'Neoplastic'

gbm = AddMetaData(gbm, metadata = df)

p = DotPlot_scCustom(gbm, features = sravya_genes, 
                     group.by = 'anno',
                     #idents = c('Myeloid', 'Neoplastic', 'Oligo', 'OPC', 'Vascular'),
                     #flip_axes = T, 
                     x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_wna_NAG_sravya.pdf', plot = p, 
       device = 'pdf', height = 4, width = 5, units = 'in')

df = data.frame(sample = gbmap@assays$RNA@data@Dimnames[[2]], 
                anno = as.character(gbmap$annotation_level_3))
df[df$anno == 'AC-like', 'anno'] = 'Neoplastic'
df[df$anno == 'NPC-like', 'anno'] = 'Neoplastic'
df[df$anno == 'OPC-like', 'anno'] = 'Neoplastic'
df[df$anno == 'MES-like', 'anno'] = 'Neoplastic'
rownames(df) = df$sample

gbmap = AddMetaData(gbmap, metadata = df)

p = DotPlot_scCustom(gbmap, features = c('ENSG00000161653', 'ENSG00000185818',
                                         'ENSG00000108381', 'ENSG00000243989',
                                         'ENSG00000154930', 'ENSG00000131069'), 
                     group.by = 'anno',
                     #idents = c('Myeloid', 'Neoplastic', 'Oligo', 'OPC', 'Vascular'),
                     #flip_axes = T, 
                     x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_gbmap_NAG_sravya.pdf', plot = p, 
       device = 'pdf', height = 5.5, width = 5, units = 'in')

abdelfattah = readRDS('../../../GSE182109_NatCommun2022/GSE182109_RAW/GSE182109_75389_CKPA.rds')

p = DotPlot_scCustom(abdelfattah, features = sravya_genes, 
                     group.by = 'Assignment',
                     #idents = c('Myeloid', 'Neoplastic', 'Oligo', 'OPC', 'Vascular'),
                     #flip_axes = T, 
                     x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_abdelfattah_NAG_sravya.pdf', plot = p, 
       device = 'pdf', height = 4, width = 5, units = 'in')


magic_obj = magic(t(GetAssayData(seuobj, slot = 'data')))
magic_res = magic_obj$result

seuobj@meta.data[['NAGS_magic']] = magic_res$NAGS
seuobj@meta.data[['ACY1_magic']] = magic_res$ACY1

p = VlnPlot(seuobj, features = c('ACY1'), 
            group.by = 'major_cell_type',
            #cols = darmanis_colors,
) + 
  NoLegend() +
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white") #+
#geom_hline(yintercept = 0, linetype = "dashed", color = 'red')

ggsave(filename = './sravya_NAG/vlnplot_major_cell_type_ACY1.pdf', plot = p, 
       device = 'pdf', height = 4, width = 4, units = 'in')
