# Visualize output of 13C-scMFA

# Import Libraries ------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(dplyr)
library(introdataviz)
library(tidyr)

# Serine model ----------------------------------------------------------------
## Visualize output of scMFA for each patient ----------------------------------
### Read single cell fluxes ----
flux_file = 'P1N_3c_scRNA_fg_20240531-152149'

flux = read.csv(file = paste('./output_serine/', 
                             flux_file, '.csv', sep = ''), 
                row.names = 1)

### Read cell IDs and cell types ----
cellnames_types = 'P1N_cellnames_types'
cellnames_types = read.table(file = paste('./patient_input_serine/patient_celltypes/', 
                                          cellnames_types, '.txt', sep = ''), 
                row.names = 1)
astrocyte_ids = cellnames_types[cellnames_types$cell_type == 'Astrocyte', 'ids']
neuron_ids = cellnames_types[cellnames_types$cell_type == 'Neuron', 'ids']
neoplastic_ids = cellnames_types[cellnames_types$cell_type == 'Neoplastic', 'ids']

### assign cell types to single cell fluxes ----
flux[['cell_type']] = NA
flux[astrocyte_ids, 'cell_type'] = 'Astrocyte'
flux[neuron_ids, 'cell_type'] = 'Neuron'
flux[neoplastic_ids, 'cell_type'] = 'Neoplastic'
flux$cell_type = factor(flux$cell_type, 
                        levels = c('Astrocyte', 'Neuron', 'Neoplastic'))

### convert a block diagonal flux matrix to a dense matrix ----
add_com_flux = function(flux, com_flux, astro_flux, glioma_flux, neuron_flux){
  flux[[com_flux]] = NA
  flux[astrocyte_ids, com_flux] = flux[astrocyte_ids, astro_flux]
  flux[neoplastic_ids, com_flux] = flux[neoplastic_ids, glioma_flux]
  flux[neuron_ids, com_flux] = flux[neuron_ids, neuron_flux]
  return(flux)
}

flux = add_com_flux(flux, 'PG_SER', 'PGa_SERa', 'PGg_SERg', 'PGn_SERn')
flux = add_com_flux(flux, 'SER_out', 'SERa_SERaout', 'SERg_SERgout', 'SERn_SERnout')
flux = add_com_flux(flux, 'SER_tran', 'SERa_SERs', 'SERs_SERg', 'SERs_SERn')
flux = add_com_flux(flux, 'SER_GLY', 'SERa_GLYa', 'SERg_GLYg', 'SERn_GLYn')
flux = add_com_flux(flux, 'SERu_SER', 'SER0_SERa', 'SER0_SERg', 'SER0_SERn')
flux = add_com_flux(flux, 'SERp_SER', 'SERp_SERa', 'SERp_SERg', 'SERp_SERn')

### violin plots of single cell fluxes of cell types ----
# compare scfluxes between cell types
my_comparisons = list( c("Astrocyte", "Neoplastic"), 
                       c("Neoplastic", "Neuron"), 
                       c("Astrocyte", "Neuron") )

compare_scfluxes = function(flux_name){
  p = ggplot(flux, aes_string(x = 'cell_type', y = flux_name, fill = 'cell_type')) +
    geom_violin(trim = F, scale = 'width') +
    geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
    stat_compare_means(comparisons = my_comparisons, 
                       label.y = c(0, 0.05, 0.1), size = 2) +
    theme_classic(base_size = 14) +
    theme(legend.position="none") +
    scale_fill_manual(values=c("#338333", "#007BBB", "#F28500"))
  
  ggsave(filename = paste(flux_file,'vlnplot', flux_name, '.pdf', 
                          sep = '_'), 
         plot = p, device = 'pdf', width = 2.5, height = 3)
}

compare_scfluxes('PG_SER')
compare_scfluxes('SER_out')
compare_scfluxes('SER_tran')
compare_scfluxes('SER_GLY')
compare_scfluxes('SERu_SER')
compare_scfluxes('SERp_SER')

### violin plots of scfluxes for neoplastic cells ----

comparisons_1 = list(c("PG_SER", "SERp_SER"), 
                     c("PG_SER", "SERu_SER"),
                     c("PG_SER", "SER_tran"),
                     c("SERp_SER", "SERu_SER"), 
                     c("SERp_SER", "SER_tran"),
                     c("SERu_SER", "SER_tran"))

plot_violin_fluxes_celltype = function(flux, flux_file, celltype, my_comparisons){
  sub_in_flux = flux[flux$cell_type == celltype, 
                     c('PG_SER', 'SERp_SER', 'SER_tran', 'SERu_SER',  
                       'SER_out', 'SER_GLY')]
  mean_neo_fluxes = apply(sub_in_flux, 2, mean)
  sd_neo_fluxes = apply(sub_in_flux, 2, sd)
  
  write.table(mean_neo_fluxes, paste(flux_file, celltype, 'mean.txt', sep = '_'))
  write.table(sd_neo_fluxes, paste(flux_file, celltype, 'sd.txt', sep = '_'))
  
  sub_in_long = tidyr::pivot_longer(sub_in_flux, cols = 1:6)
  sub_in_long[['col']] = factor(sub_in_long$name, 
                                levels = c('PG_SER', 'SERp_SER', 
                                           'SER_tran', 'SERu_SER', 
                                           'SER_out', 'SER_GLY'))
  
  p = ggplot(sub_in_long, aes_string(x = 'col', y = 'value', fill = 'col')) +
    geom_violin(trim = F, scale = 'width') +
    geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
    stat_compare_means(comparisons = my_comparisons, 
                       label.y = c(0, 0.05, 0.1, 0, 0.05, 0.1), size = 2) +
    stat_summary(fun=mean, geom="point", size=1, shape = 4, color="#FF00FF") +
    theme_classic(base_size = 14) +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    scale_fill_manual(values=c("#25783c", "#F16623", "#975aa4", "#656263",
                                        "#d52928", "#e26969"))
                                        
  ggsave(filename = paste(flux_file, 'vlnplot', celltype, 'fluxes', '.pdf', 
                          sep = '_'), 
         plot = p, device = 'pdf', width = 4, height = 3)
}

plot_violin_fluxes_celltype(flux, flux_file, 'Neoplastic', comparisons_1)

#### violin plot of main serine sources for neoplastic cells ----
sub_in_flux = flux[flux$cell_type == 'Neoplastic', 
                   c('PG_SER', 'SERp_SER', 'SER_tran')]

sub_in_long = pivot_longer(sub_in_flux, cols = 1:3)
sub_in_long[['col']] = factor(sub_in_long$name, 
                              levels = c('PG_SER', 'SERp_SER', 
                                         'SER_tran'))

my_comparisons = list(c("PG_SER", "SERp_SER"), 
                      c("PG_SER", "SER_tran"),
                      c("SERp_SER", "SER_tran"))

p = ggplot(sub_in_long, aes_string(x = 'col', y = 'value', fill = 'col')) +
  geom_violin(trim = T, scale = 'width') +
  geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, 
                     label.y = c(0, 0.05, 0.1), size = 2) +
  stat_summary(fun=mean, geom="point", size=1, shape = 4, color="#FF00FF") +
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_fill_manual(values=c("#25783c", "#F16623", "#975aa4"))


ggsave(filename = paste(flux_file, 'vlnplot', 'neoplastic_serine_sources', '.pdf', 
                        sep = '_'), 
       plot = p, device = 'pdf', width = 2.5, height = 3)

### violin plots of scfluxes for neurons ----
plot_violin_fluxes_celltype(flux, flux_file, 'Neuron', comparisons_1)

### violin plots of scfluxes for astrocytes ----
comparisons_2 = list(c("PG_SER", "SERp_SER"), 
                      c("PG_SER", "SERu_SER"),
                      c("SERp_SER", "SERu_SER"), 
                      c("SER_tran", "SER_out"),
                      c("SER_tran", "SER_GLY"),
                      c("SER_out", "SER_GLY"))
plot_violin_fluxes_celltype(flux, flux_file, 'Astrocyte', comparisons_2)


### violin plot of mean of astrocyte and neuron fluxes from two runs of N & E ----
flux_file_E = 'P7E_3c_scRNA_fg_20240510-093245'
flux_E = read.csv(file = paste('./output_serine/', 
                               flux_file_E, '.csv', sep = ''), 
                row.names = 1)

flux_file_N = 'P7N_3c_scRNA_fg_20240510-095532'
flux_N = read.csv(file = paste('./output_serine/', 
                               flux_file_N, '.csv', sep = ''), 
                  row.names = 1)

cellnames_types = 'P7E_cellnames_types' # cell types of E or N, astrocyte and neurons are the same
cellnames_types = read.table(file = paste('./patient_input_serine/patient_celltypes/', 
                                            cellnames_types, '.txt', sep = ''), 
                             row.names = 1)
astrocyte_ids = cellnames_types[cellnames_types$cell_type == 'Astrocyte', 'ids']
neuron_ids = cellnames_types[cellnames_types$cell_type == 'Neuron', 'ids']

flux_E[['cell_type']] = NA
flux_E[astrocyte_ids, 'cell_type'] = 'Astrocyte'
flux_E[neuron_ids, 'cell_type'] = 'Neuron'

flux_N[['cell_type']] = NA
flux_N[astrocyte_ids, 'cell_type'] = 'Astrocyte'
flux_N[neuron_ids, 'cell_type'] = 'Neuron'

flux_E = flux_E[!is.na(flux_E$cell_type),]
flux_N = flux_N[!is.na(flux_N$cell_type),]

combine_astro_neur_f = function(fE, fN, new_str_f, old_str_f_a, old_str_f_n){
  fE[[new_str_f]] = NA
  fE[astrocyte_ids, new_str_f] = (fE[astrocyte_ids, old_str_f_a] + fN[astrocyte_ids, old_str_f_a]) / 2
  fE[neuron_ids, new_str_f] = (fE[neuron_ids, old_str_f_n] + fN[neuron_ids, old_str_f_n]) / 2
  return(fE)
}

flux_E = combine_astro_neur_f(flux_E, flux_N, 'PG_SER', 'PGa_SERa', 'PGn_SERn')
flux_E = combine_astro_neur_f(flux_E, flux_N, 'SER_tran', 'SERa_SERs', 'SERs_SERn')
flux_E = combine_astro_neur_f(flux_E, flux_N, 'SER_out', 'SERa_SERaout', 'SERn_SERnout')
flux_E = combine_astro_neur_f(flux_E, flux_N, 'SER_GLY', 'SERa_GLYa', 'SERn_GLYn')
flux_E = combine_astro_neur_f(flux_E, flux_N, 'SERu_SER', 'SER0_SERa', 'SER0_SERn')
flux_E = combine_astro_neur_f(flux_E, flux_N, 'SERp_SER', 'SERp_SERa', 'SERp_SERn')

# fluxes for Neuron
plot_violin_fluxes_celltype(flux_E, flux_file_E, 'Neuron', comparisons_1)

# fluxes for Astrocyte
plot_violin_fluxes_celltype(flux_E, flux_file_E, 'Astrocyte', comparisons_2)


## Serine scfluxes for all patient ---------------------------------------------
patient_scMFA = c(
  'P1E_3c_scRNA_fg_20240508-164955',
  'P1N_3c_scRNA_fg_20240508-165955',
  'P2E_3c_scRNA_fg_20240508-171302',
  'P4E_3c_scRNA_fg_20240509-171827',
  'P4N_3c_scRNA_fg_20240509-230404',
  'P6E_3c_scRNA_fg_20240510-010521',
  'P6N_3c_scRNA_fg_20240510-082306',
  'P7E_3c_scRNA_fg_20240510-093245',
  'P7N_3c_scRNA_fg_20240510-095532',
  'P8N_3c_scRNA_fg_20240513-200641')

celltype_scMFA = list.files('./patient_input_serine/patient_celltypes/')
patient_ids = c('P1E', 'P1N', 'P2E', 'P4E', 'P4N', 
                'P6E', 'P6N', 'P7E', 'P7N', 'P8N')

# combine all scflux files of all patients
neo_flux = c()
for (ps in 1:length(patient_scMFA)){
  flux = read.csv(file = paste('./output_serine/', 
                               patient_scMFA[ps], '.csv', sep = ''), 
                  row.names = 1)
  
  cellnames_types = read.table(file = paste('./patient_input_serine/patient_celltypes/', 
                                            celltype_scMFA[ps], sep = ''), 
                               row.names = 1)
  
  astrocyte_ids = cellnames_types[cellnames_types$cell_type == 'Astrocyte', 'ids']
  neuron_ids = cellnames_types[cellnames_types$cell_type == 'Neuron', 'ids']
  neoplastic_ids = cellnames_types[cellnames_types$cell_type == 'Neoplastic', 'ids']
  
  flux[['cell_type']] = NA
  flux[astrocyte_ids, 'cell_type'] = 'Astrocyte'
  flux[neuron_ids, 'cell_type'] = 'Neuron'
  flux[neoplastic_ids, 'cell_type'] = 'Neoplastic'
  
  flux$cell_type = factor(flux$cell_type, levels = c('Astrocyte', 'Neuron', 'Neoplastic'))
  
  flux = add_com_flux(flux, 'PG_SER', 'PGa_SERa', 'PGg_SERg', 'PGn_SERn')
  flux = add_com_flux(flux, 'SER_out', 'SERa_SERaout', 'SERg_SERgout', 'SERn_SERnout')
  flux = add_com_flux(flux, 'SER_tran', 'SERa_SERs', 'SERs_SERg', 'SERs_SERn')
  flux = add_com_flux(flux, 'SER_GLY', 'SERa_GLYa', 'SERg_GLYg', 'SERn_GLYn')
  flux = add_com_flux(flux, 'SERu_SER', 'SER0_SERa', 'SER0_SERg', 'SER0_SERn')
  flux = add_com_flux(flux, 'SERp_SER', 'SERp_SERa', 'SERp_SERg', 'SERp_SERn')
  
  flux[['patient_id']] = rep(patient_ids[ps], dim(cellnames_types)[1]) # number of cells
  
  neo_flux = c(neo_flux, list(flux))
}

names(neo_flux) = patient_ids

neo_flux_df = bind_rows(neo_flux)

### astrocyte serine sources ----
sub_in_flux = neo_flux_df[neo_flux_df$cell_type == 'Astrocyte', 
                          c('PG_SER', 'SERp_SER','patient_id')] 

sub_in_flux[['patient_cortex']] = substr(sub_in_flux$patient_id, 1, 2)

sub_in_long = tidyr::pivot_longer(sub_in_flux, cols = 1:2) 
sub_in_long[['col']] = factor(sub_in_long$name, 
                              levels = c('PG_SER', 'SERp_SER')) 

stats_serine_sources = compare_means(data = sub_in_long, 
                                     formula = value ~ name, 
                                     group.by = 'patient_cortex', 
                                     p.adjust.method = 'holm')
write.csv(stats_serine_sources, 
          file = "stats_serine_sources_wilcox_holm_astrocyte.csv")

p = ggplot(sub_in_long, aes_string(x = 'patient_cortex', y = 'value', fill = 'col')) +
  geom_violin(trim = T, scale = 'width', position = position_dodge(0.7), 
              alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               position = position_dodge(0.7)) +
  stat_summary(fun=mean, geom="point", size=1, shape = 1, 
               position = position_dodge(0.7), color="#FF00FF") + 
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_fill_manual(values=c("#25783c", "#F16623"))  

ggsave(filename = 'all_patients_serine_sources_astrocyte_flux.pdf',
       plot = p, device = 'pdf', width = 4.5, height = 3)

### neuron serine sources ----
sub_in_flux = neo_flux_df[neo_flux_df$cell_type == 'Neuron', 
                          c('PG_SER', 'SERp_SER', 'SER_tran', 'patient_id')]

sub_in_flux[['patient_cortex']] = substr(sub_in_flux$patient_id, 1, 2)

sub_in_long = tidyr::pivot_longer(sub_in_flux, cols = 1:3) 
sub_in_long[['col']] = factor(sub_in_long$name, 
                              levels = c('PG_SER', 'SERp_SER', 'SER_tran')) 

stats_serine_sources = compare_means(data = sub_in_long, 
                                     formula = value ~ name, 
                                     group.by = 'patient_cortex', 
                                     p.adjust.method = 'holm')
write.csv(stats_serine_sources, 
          file = "stats_serine_sources_wilcox_holm_neuron.csv")

p = ggplot(sub_in_long, aes_string(x = 'patient_cortex', y = 'value', fill = 'col')) +
  geom_violin(trim = T, scale = 'width', position = position_dodge(0.7), 
              alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               position = position_dodge(0.7)) +
  stat_summary(fun=mean, geom="point", size=1, shape = 1, 
               position = position_dodge(0.7), color="#FF00FF") + 
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_fill_manual(values=c("#25783c", "#F16623", "#975aa4"))

ggsave(filename = 'all_patients_serine_sources_neuron_flux.pdf',
       plot = p, device = 'pdf', width = 5.5, height = 3)

### neoplastic serine sources ----
sub_in_flux = neo_flux_df[neo_flux_df$cell_type == 'neoplastic', 
                          c('PG_SER', 'SERp_SER', 'SER_tran', 'patient_id')]

sub_in_long = tidyr::pivot_longer(sub_in_flux, cols = 1:3) 
sub_in_long[['col']] = factor(sub_in_long$name, 
                              levels = c('PG_SER', 'SERp_SER', 'SER_tran'))

stats_serine_sources = compare_means(data = sub_in_long, 
                                     formula = value ~ name, 
                                     group.by = 'patient_id', 
                                     p.adjust.method = 'holm')
write.csv(stats_serine_sources, 
          file = "stats_serine_sources_wilcox_holm_neoplastic.csv")

p = ggplot(sub_in_long, aes_string(x = 'patient_id', y = 'value', fill = 'col')) +
  geom_violin(trim = T, scale = 'width', position = position_dodge(0.7), 
              alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               position = position_dodge(0.7)) +

  stat_summary(fun=mean, geom="point", size=1, shape = 1, 
               position = position_dodge(0.7), color="#FF00FF") + 
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_fill_manual(values=c("#25783c", "#F16623", "#975aa4"))


ggsave(filename = 'all_patients_serine_sources_neoplastic_flux.pdf',
       plot = p, device = 'pdf', width = 9, height = 3)

### serine consumption across cell types ----
sub_in_long = neo_flux_df[, 
                          c('SER_GLY' , 'SER_out', 'cell_type',
                            'patient_id')]

sub_in_long[sub_in_long$SER_GLY < 0, 'SER_GLY'] = 0
sub_in_long[sub_in_long$SER_out < 0, 'SER_out'] = 0

sub_in_long[['SER_GLY_SER_cons']] = sub_in_long$SER_GLY / (sub_in_long$SER_GLY + 
                                                             sub_in_long$SER_out)

sub_in_long[['SER_GLY_SER_cons']] = sub_in_long$SER_GLY_SER_cons * 100
sub_in_long[is.na(sub_in_long$SER_GLY_SER_cons), 'SER_GLY_SER_cons'] = 0
apply(sub_in_long, 2, max)

sub_in_long[['SER_GLY']] = NULL
sub_in_long[['SER_out']] = NULL

sub_in_long[['col']] = factor(sub_in_long$cell_type, 
                              levels = c('Astrocyte', 'Neuron', 'Neoplastic'))


my_comparisons = list(c('Neuron', 'Noeplastic'),
                      c('Astrocyte', 'Neoplastic'))
p = ggplot(sub_in_long, aes_string(x = 'patient_id', y = 'SER_GLY_SER_cons', 
                                   fill = 'col')) +
  geom_violin(trim = T, scale = 'width', position = position_dodge(0.7), 
              alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               position = position_dodge(0.7)) +
  stat_summary(fun=mean, geom="point", size=1, shape = 1, 
               position = position_dodge(0.7), color="#FF00FF") + 
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_fill_manual(values=c("#1a741a", "#3395c9", "#f28500")) 

ggsave(filename = 'all_patients_ser_gly_ser_cons_flux_celltype_stat.pdf',
       plot = p, device = 'pdf', width = 9, height = 3)

stats_ser_gly_ser_cons = compare_means(data = sub_in_long, 
                                       formula = SER_GLY_SER_cons ~ cell_type, 
                                       group.by = 'patient_id', 
                                       ref.group = 'Neoplastic', 
                                       p.adjust.method = 'holm')

write.csv(stats_ser_gly_ser_cons, file = "stats_ser_gly_ser_cons_wilcox_holm.csv")

## Correlation of 13C-scMFA estimated fluxes with CNN predicted fluxes ---------
scMFA_dir = "./output_serine"
patient_scMFA = c(
  'P1E_3c_scRNA_fg_20240508-164955',
  'P1N_3c_scRNA_fg_20240508-165955',
  'P2E_3c_scRNA_fg_20240508-171302',
  'P4E_3c_scRNA_fg_20240509-171827',
  'P4N_3c_scRNA_fg_20240509-230404',
  'P6E_3c_scRNA_fg_20240510-010521',
  'P6N_3c_scRNA_fg_20240510-082306',
  'P7E_3c_scRNA_fg_20240510-093245',
  'P7N_3c_scRNA_fg_20240510-095532',
  'P8N_3c_scRNA_fg_20240513-200641')

neo_flux = c() 
for (ps in 1:length(patient_scMFA)){
  # read average scMFA fluxes of neoplastic cells saved & calculated in 
  # previous sections
  mal_mean = read.delim(file = paste(scMFA_dir, '/',
                                     patient_scMFA[ps], 
                                     '_Neoplastic_mean.txt', 
                                     sep = ''), sep = ' ',
                        quote = '', row.names = 1, header = T)  
  neo_flux = c(neo_flux, list(mal_mean[c(1, 2, 3), ]))
}

neo_flux_df = data.frame(neo_flux)
colnames(neo_flux_df) = c('P1E', 'P1N', 'P2E', 'P4E', 'P4N', 
                          'P6E', 'P6N', 'P7E', 'P7N', 'P8N')
rownames(neo_flux_df) = c('denovo', 'plasma', 'tme')

neo_flux_normalized = sweep(neo_flux_df, 2, colSums(neo_flux_df),`/`)
neo_flux_normalized = data.frame(t(neo_flux_normalized))

# Directory of CNN predicted serine fluxes in patients
ML_pred_dir = "../metabolic_CNN/combined_patient_pred_serine_CNN/"

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

# remove patients who don't have scRNA-seq data
glioma_denovo[, c('Patient2N', 'Patient3E',
                  'Patient3N', 'Patient5E',
                  'Patient5N')] = list(NULL)
glioma_plasma[, c('Patient2N', 'Patient3E',
                  'Patient3N', 'Patient5E',
                  'Patient5N')] = list(NULL)
glioma_tme[, c('Patient2N', 'Patient3E',
               'Patient3N', 'Patient5E',
               'Patient5N')] = list(NULL)

colnames(glioma_denovo) = colnames(neo_flux_df)

colnames(glioma_plasma) = colnames(neo_flux_df)

colnames(glioma_tme) = colnames(neo_flux_df)

glioma_denovo_ps = colMeans(glioma_denovo)
glioma_plasma_ps = colMeans(glioma_plasma)
glioma_tme_ps = colMeans(glioma_tme)

glioma_ps = rbind(glioma_denovo_ps, glioma_plasma_ps, glioma_tme_ps)

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

# scatter plots with regression line
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


# Purine model ----
## Visualize output of scMFA for each patient ---------------------------------
# list of scMFA output patient files
patient_flux_files = c(
  'P1E_purine_scRNA_fg_20240704-154709',
  'P1N_purine_scRNA_fg_20240704-165203',
  'P2E_purine_scRNA_fg_20240704-180550',
  'P4E_purine_scRNA_fg_20240705-192842',
  'P4N_purine_scRNA_fg_20240707-083113',
  'P6E_purine_scRNA_fg_20240707-164950',
  'P6N_purine_scRNA_fg_20240707-221203',
  'P7E_purine_scRNA_fg_20240708-001347',
  'P7N_purine_scRNA_fg_20240708-014219',
  'P8N_purine_scRNA_fg_20240708-212653')

patient_ids = substr(patient_flux_files, 1, 3)
patient_sites = substr(patient_flux_files, 3, 3)

combine_myeloid_neo_f = function(f, new_str_f, old_str_f_m, old_str_f_n){
  f[[new_str_f]] = NA
  f[myeloid_ids, new_str_f] = f[myeloid_ids, old_str_f_m]
  f[neoplastic_ids, new_str_f] = f[neoplastic_ids, old_str_f_n]
  return(f)
}

# Plot scfluxes for each patient
for (ps in 1:length(patient_flux_files)){
  flux_file = patient_flux_files[ps]
  flux = read.csv(file = paste('./output_purine/', flux_file, '.csv', sep = ''), 
                  row.names = 1)
  
  cellnames_types = paste(patient_ids[ps], '_purine_cellnames_types', sep = '')
  
  cellnames_types = read.table(file = paste('./patient_input_purine/patient_celltypes/', 
                                            cellnames_types, '.txt', sep = ''), 
                               row.names = 1)
  
  myeloid_ids = cellnames_types[cellnames_types$cell_type == 'Myeloid', 'ids']
  neoplastic_ids = cellnames_types[cellnames_types$cell_type == 'Neoplastic', 'ids']
  
  flux = combine_myeloid_neo_f(flux, 'denovo_IMP', 'denovo_IMP_m', 'denovo_IMP_g')
  flux = combine_myeloid_neo_f(flux, 'HPX_IMP', 'HPX_IMP_m', 'HPX_IMP_g')
  flux = combine_myeloid_neo_f(flux, 'IMP_INO', 'IMP_INO_m', 'IMP_INO_g')
  flux = combine_myeloid_neo_f(flux, 'IMP_AMP', 'IMP_AMP_m', 'IMP_AMP_g')
  flux = combine_myeloid_neo_f(flux, 'AMP_IMP', 'AMP_IMP_m', 'AMP_IMP_g')
  flux = combine_myeloid_neo_f(flux, 'IMP_GMP', 'IMP_GMP_m', 'IMP_GMP_g')
  flux = combine_myeloid_neo_f(flux, 'ADE_AMP', 'ADE_AMP_m', 'ADE_AMP_g')
  flux = combine_myeloid_neo_f(flux, 'AMP_out', 'AMP_out_m', 'AMP_out_g')
  flux = combine_myeloid_neo_f(flux, 'ADE_INO', 'ADE_INO_m', 'ADE_INO_g')
  flux = combine_myeloid_neo_f(flux, 'INO_HPX', 'INO_HPX_m', 'INO_HPX_g')
  flux = combine_myeloid_neo_f(flux, 'HPX_INO', 'HPX_INO_m', 'HPX_INO_g')
  flux = combine_myeloid_neo_f(flux, 'GMP_GDP', 'GMP_GDP_m', 'GMP_GDP_g')
  flux = combine_myeloid_neo_f(flux, 'GMP_GUO', 'GMP_GUO_m', 'GMP_GUO_g')
  flux = combine_myeloid_neo_f(flux, 'GUA_GMP', 'GUA_GMP_m', 'GUA_GMP_g')
  flux = combine_myeloid_neo_f(flux, 'GUO_GUA', 'GUO_GUA_m', 'GUO_GUA_g')
  flux = combine_myeloid_neo_f(flux, 'GUA_GUO', 'GUA_GUO_m', 'GUA_GUO_g')
  flux = combine_myeloid_neo_f(flux, 'GDP_out', 'GDP_out_m', 'GDP_out_g')
  
  sub_in_flux = flux[, c(
    'denovo_IMP',
    'HPX_IMP',
    'IMP_INO',
    'IMP_AMP',
    'AMP_IMP',
    'IMP_GMP',
    'ADE_AMP',
    'AMP_out',
    'ADE_INO',
    'INO_HPX',
    'HPX_INO',
    'GMP_GDP',
    'GMP_GUO',
    'GUA_GMP',
    'GUO_GUA',
    'GUA_GUO',
    'GDP_out')]
  
  # replace negative fluxes with NA
  sub_in_flux = replace(sub_in_flux, sub_in_flux < 0, NA)
  # remove cells with negative fluxes
  #sub_in_flux = subset(sub_in_flux,(rowSums(sign(sub_in_flux) < 0) == 0) )
  
  sub_in_flux[['cell_type']] = NA
  sub_in_flux[rownames(sub_in_flux) %in% myeloid_ids, 'cell_type'] = 'Myeloid'
  sub_in_flux[rownames(sub_in_flux) %in% neoplastic_ids, 'cell_type'] = 'Neoplastic'
  sub_in_flux$cell_type = factor(sub_in_flux$cell_type, 
                                 levels = c('Myeloid', 'Neoplastic'))
  
  sub_in_long = tidyr::pivot_longer(sub_in_flux, cols = 1:17)
  sub_in_long$name = factor(sub_in_long$name, 
                            levels = c('denovo_IMP',
                                       'HPX_IMP',
                                       'IMP_INO',
                                       'IMP_AMP',
                                       'AMP_IMP',
                                       'IMP_GMP',
                                       'ADE_AMP',
                                       'AMP_out',
                                       'ADE_INO',
                                       'INO_HPX',
                                       'HPX_INO',
                                       'GMP_GDP',
                                       'GMP_GUO',
                                       'GUA_GMP',
                                       'GUO_GUA',
                                       'GUA_GUO',
                                       'GDP_out'))
  
  p = ggplot(sub_in_long, aes(x = name, y = value, fill = cell_type)) +
    introdataviz::geom_split_violin(alpha = 0.7, scale = 'width') +
    geom_boxplot(outlier.shape = NA, width = 0.3,
                 show.legend = FALSE) +
    scale_x_discrete(name = "Purine Fluxes") +
    scale_fill_manual(values = c("#e7f582", "#f28500"))+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90))
  
  ggsave(filename = paste(flux_file,
                          'vlnplot_myeloid_neoplastic_purine_fluxes.pdf',
                          sep = '_'),
         plot = p, device = 'pdf', height = 4, width = 8, units = 'in')
  
  stats_gmp_sources = compare_means(data = sub_in_long,
                                    formula = value ~ name,
                                    group.by = 'cell_type',
                                    p.adjust.method = 'holm')
  write.csv(stats_gmp_sources,
            file = paste(flux_file,
                         "_stats_gmp_sources_wilcox_holm.csv", sep = ''))
  
  stats_gmp_sources = compare_means(data = sub_in_long, 
                                    formula = value ~ cell_type, 
                                    group.by = 'name', 
                                    p.adjust.method = 'holm')
  write.csv(stats_gmp_sources, 
            file = paste(flux_file, 
                         "_stats_purine_compare_neo_mye_wilcox_holm.csv", 
                         sep = ''))
  
  flux_myeloid = sub_in_flux[sub_in_flux$cell_type == 'Myeloid',
                             !colnames(sub_in_flux) %in% c('cell_type')]
  flux_neoplastic = sub_in_flux[sub_in_flux$cell_type == 'Neoplastic',
                                !colnames(sub_in_flux) %in% c('cell_type')]
  
  flux_myeloid = na.omit(flux_myeloid)
  flux_neoplastic = na.omit(flux_neoplastic)
  
  mean_myeloid_fluxes = apply(flux_myeloid, 2, mean)
  mean_neoplastic_fluxes = apply(flux_neoplastic, 2, mean)
  
  write.table(mean_myeloid_fluxes, paste(flux_file, 'myeloid_mean.txt'))
  write.table(mean_neoplastic_fluxes, paste(flux_file, 'neoplastic_mean.txt'))
}  


## Purine scfluxes for all patients  -------------------------------------------
patient_flux_files = c(
  'P1E_purine_scRNA_fg_20240704-154709',
  'P1N_purine_scRNA_fg_20240704-165203',
  'P2E_purine_scRNA_fg_20240704-180550',
  'P4E_purine_scRNA_fg_20240705-192842',
  'P4N_purine_scRNA_fg_20240707-083113',
  'P6E_purine_scRNA_fg_20240707-164950',
  'P6N_purine_scRNA_fg_20240707-221203',
  'P7E_purine_scRNA_fg_20240708-001347',
  'P7N_purine_scRNA_fg_20240708-014219',
  'P8N_purine_scRNA_fg_20240708-212653')

celltype_scMFA = list.files('./patient_input_purine/patient_celltypes/')

patient_ids = c('P1E', 'P1N', 'P2E', 'P4E', 'P4N', 
                'P6E', 'P6N', 'P7E', 'P7N', 'P8N')
neo_flux = c()
for (ps in 1:length(patient_flux_files)){
  flux = read.csv(file = paste('./output_purine/', 
                               patient_flux_files[ps], '.csv', sep = ''), 
                  row.names = 1)
  
  cellnames_types = read.table(file = paste('./patient_input_purine/patient_celltypes/', 
                                            celltype_scMFA[ps], sep = ''), 
                               row.names = 1)
  
  myeloid_ids = cellnames_types[cellnames_types$cell_type == 'Myeloid', 'ids']
  neoplastic_ids = cellnames_types[cellnames_types$cell_type == 'Neoplastic', 'ids']
  
  flux = combine_myeloid_neo_f(flux, 'denovo_IMP', 'denovo_IMP_m', 'denovo_IMP_g')
  flux = combine_myeloid_neo_f(flux, 'HPX_IMP', 'HPX_IMP_m', 'HPX_IMP_g')
  flux = combine_myeloid_neo_f(flux, 'IMP_INO', 'IMP_INO_m', 'IMP_INO_g')
  flux = combine_myeloid_neo_f(flux, 'IMP_AMP', 'IMP_AMP_m', 'IMP_AMP_g')
  flux = combine_myeloid_neo_f(flux, 'AMP_IMP', 'AMP_IMP_m', 'AMP_IMP_g')
  flux = combine_myeloid_neo_f(flux, 'IMP_GMP', 'IMP_GMP_m', 'IMP_GMP_g')
  flux = combine_myeloid_neo_f(flux, 'ADE_AMP', 'ADE_AMP_m', 'ADE_AMP_g')
  flux = combine_myeloid_neo_f(flux, 'AMP_out', 'AMP_out_m', 'AMP_out_g')
  flux = combine_myeloid_neo_f(flux, 'ADE_INO', 'ADE_INO_m', 'ADE_INO_g')
  flux = combine_myeloid_neo_f(flux, 'INO_HPX', 'INO_HPX_m', 'INO_HPX_g')
  flux = combine_myeloid_neo_f(flux, 'HPX_INO', 'HPX_INO_m', 'HPX_INO_g')
  flux = combine_myeloid_neo_f(flux, 'GMP_GDP', 'GMP_GDP_m', 'GMP_GDP_g')
  flux = combine_myeloid_neo_f(flux, 'GMP_GUO', 'GMP_GUO_m', 'GMP_GUO_g')
  flux = combine_myeloid_neo_f(flux, 'GUA_GMP', 'GUA_GMP_m', 'GUA_GMP_g')
  flux = combine_myeloid_neo_f(flux, 'GUO_GUA', 'GUO_GUA_m', 'GUO_GUA_g')
  flux = combine_myeloid_neo_f(flux, 'GUA_GUO', 'GUA_GUO_m', 'GUA_GUO_g')
  flux = combine_myeloid_neo_f(flux, 'GDP_out', 'GDP_out_m', 'GDP_out_g')
  
  sub_in_flux = flux[, c(
    'denovo_IMP',
    'HPX_IMP',
    'IMP_INO',
    'IMP_AMP',
    'AMP_IMP',
    'IMP_GMP',
    'ADE_AMP',
    'AMP_out',
    'ADE_INO',
    'INO_HPX',
    'HPX_INO',
    'GMP_GDP',
    'GMP_GUO',
    'GUA_GMP',
    'GUO_GUA',
    'GUA_GUO',
    'GDP_out')]
  
  # replace negative fluxes with NA
  sub_in_flux = replace(sub_in_flux, sub_in_flux < 0, NA)
  # remove cells with negative fluxes
  #sub_in_flux = subset(sub_in_flux,(rowSums(sign(sub_in_flux) < 0) == 0) )
  
  sub_in_flux[['cell_type']] = NA
  sub_in_flux[rownames(sub_in_flux) %in% myeloid_ids, 'cell_type'] = 'Myeloid'
  sub_in_flux[rownames(sub_in_flux) %in% neoplastic_ids, 'cell_type'] = 'Neoplastic'
  sub_in_flux$cell_type = factor(sub_in_flux$cell_type, 
                                 levels = c('Myeloid', 'Neoplastic'))
  
  sub_in_flux[['patient_id']] = rep(patient_ids[ps], dim(cellnames_types)[1]) # number of cells
  
  neo_flux = c(neo_flux, list(sub_in_flux))
}

names(neo_flux) = patient_ids

neo_flux_df = bind_rows(neo_flux)

neo_flux_df[['ratio']] = neo_flux_df$IMP_GMP / (neo_flux_df$IMP_GMP + neo_flux_df$GUA_GMP)
neo_flux_df[['GMP_synthesis']] = neo_flux_df$IMP_GMP + neo_flux_df$GUA_GMP
neo_flux_df[['net_IMP_AMP']] = neo_flux_df$IMP_AMP - neo_flux_df$AMP_IMP
neo_flux_df[['IMP_synthesis']] = neo_flux_df$denovo_IMP + neo_flux_df$IMP_GMP
neo_flux_df[['net_INO_HPX']] = neo_flux_df$INO_HPX - neo_flux_df$HPX_INO
neo_flux_df[['net_GUO_GUA']] = neo_flux_df$GUO_GUA - neo_flux_df$GUA_GUO

p = ggplot(neo_flux_df, aes(x = patient_id, y = net_GUO_GUA, fill = cell_type)) +
  introdataviz::geom_split_violin(alpha = 0.7, scale = 'width') +
  geom_boxplot(outlier.shape = NA, width = 0.3,
               show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'red') +
  scale_x_discrete(name = "Patient IDs") +
  scale_fill_manual(values = c("#e7f582", "#f28500"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.1))

ggsave(filename = paste('vlnplot_myeloid_neoplastic_net_GUO_GUA.pdf',
                        sep = '_'),
       plot = p, device = 'pdf', height = 4, width = 5.5, units = 'in')

