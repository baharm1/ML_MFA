# Visualize output of modified scFEA
# scRNA-seq dataset from Darmanis et al., Cell Reports 2017

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

# Read Darmanis dataset -------------------------------------------------------
# Load preprocessed data saved in scRNA_analysis_Darmanis.R
seuobj = readRDS('darmanis.rds')
seuobj = SetIdent(seuobj, value = c('major_cell_type'))
metadata = seuobj@meta.data

# Visualize output of modified scFEA ------------------------------------------
## Cell ids of each cell type ----
astrocyte_ids = rownames(metadata[metadata$major_cell_type == 'Astrocyte', ])
neuron_ids = rownames(metadata[metadata$major_cell_type == 'Neuron', ])
neoplastic_ids = rownames(metadata[metadata$major_cell_type == 'Neoplastic', ])

ids = c(astrocyte_ids, neuron_ids, neoplastic_ids)

## Read output of modified scFEA ----
flux_file = 'three_comp_scRNA_f_20240415-175335'
flux = read.csv(file = paste('./scfea/myEditsFor_scFEA/output/', 
                             flux_file, '.csv', sep = ''), 
                row.names = 1)

## assign cell types to single cell fluxes ----
flux[['cell_type']] = NA
flux[astrocyte_ids, 'cell_type'] = 'Astrocyte'
flux[neuron_ids, 'cell_type'] = 'Neuron'
flux[neoplastic_ids, 'cell_type'] = 'Neoplastic'
flux$cell_type = factor(flux$cell_type, 
                        levels = c('Astrocyte', 'Neuron', 'Neoplastic'))

## convert a block diagonal flux matrix to a dense matrix ----
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

## violin plots of single cell fluxes ----
# compare scfluxes between cell types
my_comparisons = list( c("Astrocyte", "Neoplastic"), 
                       c("Neoplastic", "Neuron"), 
                       c("Astrocyte", "Neuron") )

compare_scfluxes = function(flux_name){
  p = ggplot(flux, aes_string(x = 'cell_type', y = flux_name, fill = 'cell_type')) +
    geom_violin(trim = F, scale = 'width') +
    geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
    stat_compare_means(comparisons = my_comparisons, 
                       label.y = c(0.1, 0.05, 0), size = 2) +
    theme_classic(base_size = 14) +
    theme(legend.position="none") +
    scale_fill_manual(values=c("#338333", "#007BBB", "#F28500"))
  
  ggsave(filename = paste('vlnplot', flux_file, flux_name, '.pdf', 
                          sep = '_'), 
         plot = p, device = 'pdf', width = 2.5, height = 3)
}

compare_scfluxes('PG_SER')
compare_scfluxes('SER_out')
compare_scfluxes('SER_tran')
compare_scfluxes('SER_GLY')

