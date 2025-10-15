# Visualize output of 13C-scMFA

# Import Libraries ------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(dplyr)
library(introdataviz)
library(tidyr)

# Serine model validation -----------------------------------------------------
## GBM12, GBM38, HF2303 -------------------------------------------------------
mouse_scMFA = c(
  'GBM12_4c_fg_20250418-144718',
  'GBM38_4c_fg_20250418-183546',
  'HF2303_4c_fg_20250419-000940'
)

celltype_scMFA = list.files('./GBM12_GBM38_HF2303_input_serine/celltypes/')
mouse_ids = c('GBM12', 'GBM38', 'HF2303')

### convert a block diagonal flux matrix to a dense matrix ----
add_com_flux = function(flux, com_flux, astro_flux, neuron_flux, neo_flux){
  flux[[com_flux]] = NA
  flux[astrocyte_ids, com_flux] = flux[astrocyte_ids, astro_flux]
  flux[neuron_ids, com_flux] = flux[neuron_ids, neuron_flux]
  flux[inv_neoplastic_ids, com_flux] = flux[inv_neoplastic_ids, neo_flux]
  flux[neoplastic_ids, com_flux] = flux[neoplastic_ids, neo_flux]
  return(flux)
}

# combine all scflux files of all mouse models
neo_flux = c()
for (ps in 1:length(mouse_scMFA)){
  flux = read.csv(file = paste('./output_serine_mouse/', 
                               mouse_scMFA[ps], '.csv', sep = ''), 
                  row.names = 1)
  
  cellnames_types = read.table(file = paste('./GBM12_GBM38_HF2303_input_serine/celltypes/', 
                                            celltype_scMFA[ps], sep = ''), 
                               row.names = 1)
  
  astrocyte_ids = cellnames_types[cellnames_types$cell_type == 'Astrocyte', 'ids']
  neuron_ids = cellnames_types[(cellnames_types$cell_type == 'Neuron_1') |
                                 (cellnames_types$cell_type == 'Neuron_2'), 'ids']
  inv_neoplastic_ids = cellnames_types[cellnames_types$cell_type == 'Invading_Neoplastic', 'ids']
  neoplastic_ids = cellnames_types[cellnames_types$cell_type == 'Neoplastic', 'ids']
  
  flux[['cell_type']] = NA
  flux[astrocyte_ids, 'cell_type'] = 'Astrocyte'
  flux[neuron_ids, 'cell_type'] = 'Neuron'
  flux[inv_neoplastic_ids, 'cell_type'] = 'Invading_Neoplastic'
  flux[neoplastic_ids, 'cell_type'] = 'Neoplastic'
  
  flux$cell_type = factor(flux$cell_type, levels = c('Astrocyte', 'Neuron', 
                                                     'Invading_Neoplastic',
                                                     'Neoplastic'))
  
  flux = add_com_flux(flux, 'PG_SER', 'PGa_SERa', 'PGn_SERn', 'PGg_SERg')
  flux = add_com_flux(flux, 'SER_out', 'SERa_SERaout', 'SERn_SERnout', 
                      'SERg_SERgout')
  flux = add_com_flux(flux, 'SER_tran', 'SERa_SERs', 'SERs_SERn', 'SERs_SERg')
  flux = add_com_flux(flux, 'SER_GLY', 'SERa_GLYa', 'SERn_GLYn', 'SERg_GLYg')
  flux = add_com_flux(flux, 'SERu_SER', 'SER0_SERa', 'SER0_SERn', 'SER0_SERg') 
  flux = add_com_flux(flux, 'SERp_SER', 'SERp_SERa', 'SERp_SERn', 'SERp_SERg')
  
  flux[['mouse_id']] = rep(mouse_ids[ps], dim(cellnames_types)[1]) # number of cells
  
  neo_flux = c(neo_flux, list(flux))
}

names(neo_flux) = mouse_ids

neo_flux_df = bind_rows(neo_flux)

### ratio of glycine formation to total serine consumption across cell types ----
sub_in_long = neo_flux_df[, 
                          c('SER_GLY' , 'SER_out', 'cell_type',
                            'mouse_id')]

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
                              levels = c('Astrocyte', 'Neuron', 
                                         'Invading_Neoplastic', 'Neoplastic'))
# Figure 6M
p = ggplot(sub_in_long, aes_string(x = 'mouse_id', y = 'SER_GLY_SER_cons', 
                                   fill = 'col')) +
  geom_violin(trim = T, scale = 'width', position = position_dodge(0.7), 
              alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               position = position_dodge(0.7)) +
  stat_summary(fun=mean, geom="point", size=1, shape = 1, 
               position = position_dodge(0.7), color="#FF00FF") + 
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_fill_manual(values=c("#1a741a", "#3395c9", "#d9aea3", "#f28500")) 

ggsave(filename = 'all_mice_ser_gly_ser_cons_flux_celltype_stat_v3.pdf',
       plot = p, device = 'pdf', width = 5.5, height = 3)

stats_ser_gly_ser_cons = compare_means(data = sub_in_long, 
                                       formula = SER_GLY_SER_cons ~ cell_type, 
                                       group.by = 'mouse_id', 
                                       ref.group = 'Neoplastic', 
                                       p.adjust.method = 'holm')

write.csv(stats_ser_gly_ser_cons, file = "stats_ser_gly_ser_cons_wilcox_holm_v3.csv")

### serine sources across cell types ----
# Figure 6J-L
# replace 'SER_tran' with 'SERp_SER' and 'PG_SER'
sub_in_long = neo_flux_df[, 
                          c('SER_tran' , 'cell_type',
                            'mouse_id')] 

sub_in_long[['col']] = factor(sub_in_long$cell_type, 
                              levels = c('Astrocyte', 'Neuron', 
                                         'Invading_Neoplastic', 'Neoplastic'))

p = ggplot(sub_in_long, aes_string(x = 'mouse_id', y = 'SER_tran', 
                                   fill = 'col')) +
  geom_violin(trim = F, scale = 'width', position = position_dodge(0.7), 
              alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               position = position_dodge(0.7)) +
  stat_summary(fun=mean, geom="point", size=1, shape = 1, 
               position = position_dodge(0.7), color="#FF00FF") + 
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_fill_manual(values=c("#1a741a", "#3395c9", "#d9aea3", "#f28500")) 

ggsave(filename = 'all_mice_SER_trans_flux_celltype_stat_v3.pdf',
       plot = p, device = 'pdf', width = 5.5, height = 3)

stats_ser = compare_means(data = sub_in_long, 
                          formula = SER_tran ~ cell_type, 
                          group.by = 'mouse_id', 
                          ref.group = 'Neoplastic', 
                          p.adjust.method = 'holm')

write.csv(stats_ser, file = "stats_SER_trans_wilcox_holm_v3.csv")

# Purine model validation -----------------------------------------------------
## TRP vs GBM38 ---------------------------------------------------------------

TRP_GBM38_files = c('GBM38_SP1_purine_scRNA_fg_20250729-173516',
                    'TRP_SP3_purine_scRNA_fg_20250729-154613')

celltype_scMFA = c('GBM38_SP1_purine_cellnames_types.txt',
                   'TRP_SP3_purine_cellnames_types.txt')

model_ids = c('GBM38', 'TRP')

combine_myeloid_neo_f = function(f, new_str_f, old_str_f_m, old_str_f_n){
  f[[new_str_f]] = NA
  f[myeloid_ids, new_str_f] = f[myeloid_ids, old_str_f_m]
  f[neoplastic_ids, new_str_f] = f[neoplastic_ids, old_str_f_n]
  return(f)
}

neo_flux = c()
for (ps in 1:length(TRP_GBM38_files)){
  flux = read.csv(file = paste('./output_purine/', 
                               TRP_GBM38_files[ps], '.csv', sep = ''), 
                  row.names = 1)
  
  cellnames_types = read.table(file = paste('./TRP_GBM38_input_purine/celltypes/', 
                                            celltype_scMFA[ps], sep = ''), 
                               row.names = 1)
  
  myeloid_ids = cellnames_types[cellnames_types$cell_type == 'Myeloid', 'ids']
  neoplastic_ids = cellnames_types[cellnames_types$cell_type == 'Neoplastic', 'ids']
  
  flux = combine_myeloid_neo_f(flux, 'IMP_GMP', 'IMP_GMP_m', 'IMP_GMP_g')
  flux = combine_myeloid_neo_f(flux, 'GUA_GMP', 'GUA_GMP_m', 'GUA_GMP_g')
  
  sub_in_flux = flux[, c('IMP_GMP', 'GUA_GMP')]
  
  # replace negative fluxes with NA
  sub_in_flux = replace(sub_in_flux, sub_in_flux < 0, NA)
  
  sub_in_flux[['cell_type']] = NA
  sub_in_flux[rownames(sub_in_flux) %in% myeloid_ids, 'cell_type'] = 'Myeloid'
  sub_in_flux[rownames(sub_in_flux) %in% neoplastic_ids, 'cell_type'] = 'Neoplastic'
  sub_in_flux$cell_type = factor(sub_in_flux$cell_type, 
                                 levels = c('Myeloid', 'Neoplastic'))
  
  sub_in_flux[['model_id']] = rep(model_ids[ps], dim(cellnames_types)[1]) # number of cells
  
  neo_flux = c(neo_flux, list(sub_in_flux))
}

names(neo_flux) = model_ids

neo_flux_df = bind_rows(neo_flux)
neo_flux_df = neo_flux_df[!is.na(neo_flux_df$IMP_GMP), ]

# Figure 7E
p = ggplot(neo_flux_df, aes(x = model_id, y = IMP_GMP, fill = cell_type)) +
  introdataviz::geom_split_violin(alpha = 0.7, scale = 'width') +
  geom_boxplot(outlier.shape = NA, width = 0.3,
               show.legend = FALSE) +
  scale_x_discrete(name = "Mouse IDs") +
  scale_fill_manual(values = c("#e7f582", "#f28500"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.1))

ggsave(filename = 'GBM38_TRP_IMP_GMP.pdf',
       plot = p, device = 'pdf', height = 4, width = 3, units = 'in')


# stat calculation
sub_in_flux = neo_flux_df[neo_flux_df$cell_type == 'Neoplastic', 
                          c('IMP_GMP', 'ratio','model_id')] 

sub_in_long = tidyr::pivot_longer(sub_in_flux, cols = 1:2) 

stats_purine_sources = ggpubr::compare_means(data = sub_in_long, 
                                     formula = value ~ model_id, 
                                     group.by = 'name', 
                                     p.adjust.method = 'holm')
write.csv(stats_purine_sources, 
          file = "GBM38_TRP_stats_purine_sources_wilcox_holm_neoplastic.csv")
