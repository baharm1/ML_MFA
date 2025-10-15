# Prepare input data for modified scFEA
# scRNA-seq dataset from Darmanis et al., Cell Reports 2017

# Import Libraries ------------------------------------------------------------
# R 4.2.2
library(Seurat) # version 4.2.0
library(scCustomize) # version 1.0.0
library(RColorBrewer) # version 1.1-3
library(ggplot2) # version 3.4.2

# Read Darmanis dataset -------------------------------------------------------
# Load preprocessed data saved in scRNA_analysis_Darmanis.R
seuobj = readRDS('darmanis.rds')
seuobj = SetIdent(seuobj, value = c('major_cell_type'))

# Define serine related genes -------------------------------------------------
list_serine_genes = c('PHGDH', 'PSAT1', 'PSPH', 
                      'SHMT1', 'SHMT2', 
                      'SDS', 'SDSL',
                      'SRR',
                      'CBS', 
                      'SPTLC1', 'SPTLC2', 'SPTLC3', 
                      'SARS', 'SARS2', 
                      'PTDSS1', 'PTDSS2',
                      'SLC1A4', 'SLC1A5', 'SLC3A2', 
                      'SLC6A9', 'SLC6A14', 'SLC7A8', 'SLC7A10', 'SLC12A4', 
                      'SLC25A1', 'SLC25A15', 'SLC25A28',
                      'SLC36A1',
                      'SLC38A1', 'SLC38A2', 'SLC38A4', 'SLC38A5', 'SLC38A7')

# Figure S2I
p = DotPlot_scCustom(seuobj, features = list_serine_genes, 
                     group.by = 'major_cell_type',
                     flip_axes = T, x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_major_cell_type_serine_all.pdf', plot = p, 
       device = 'pdf', height = 10, width = 6, units = 'in')
	   
# Create block diagonal matrix ------------------------------------------------
# normalized gene expression of each cell type are combined in 
# a block diagonal matrix: 
# O - -
# - O -
# - - O
# Where O is the normalized gene expression matrix of a cell type with 
# the shape of number of genes by number of cells of that cell type

astro = subset(seuobj, idents = c('Astrocyte'))
astro_ge = GetAssayData(astro, slot = 'data')
astro_ge = astro_ge[list_serine_genes, ]
dim(astro_ge) # 31 * 88
dimnames(astro_ge)[[1]] = paste(dimnames(astro_ge)[[1]], 'a', sep = '_')
astro_ge = as.matrix(astro_ge)

neuron = subset(seuobj, idents = c('Neuron'))
neuron_ge = GetAssayData(neuron, slot = 'data')
neuron_ge = neuron_ge[list_serine_genes, ]
dim(neuron_ge) # 31 * 21
dimnames(neuron_ge)[[1]] = paste(dimnames(neuron_ge)[[1]], 'n', sep = '_')
neuron_ge = as.matrix(neuron_ge)

neoplastic = subset(seuobj, idents = c('Neoplastic'))
neoplastic_ge = GetAssayData(neoplastic, slot = 'data')
neoplastic_ge = neoplastic_ge[list_serine_genes, ]
dim(neoplastic_ge) # 31 * 1091
dimnames(neoplastic_ge)[[1]] = paste(dimnames(neoplastic_ge)[[1]], 'g', sep = '_')
neoplastic_ge = as.matrix(neoplastic_ge)

astro_ge = as.data.frame(astro_ge)
neuron_ge = as.data.frame(neuron_ge)
neoplastic_ge = as.data.frame(neoplastic_ge)

combined = bind_rows(astro_ge, neuron_ge, neoplastic_ge)
combined[is.na(combined)] = 0

write.table(combined, file = 'three_comp_scRNA.txt', sep = '\t', quote = F)
write.csv(combined, file = 'three_comp_scRNA.csv', quote = F)

cellnames_types = data.frame(ids, row.names = ids)
cellnames_types[['cell_type']] = NA
cellnames_types[cellnames_types$ids %in% astrocyte_ids, 'cell_type'] = 'Astrocyte'
cellnames_types[cellnames_types$ids %in% neuron_ids, 'cell_type'] = 'Neuron'
cellnames_types[cellnames_types$ids %in% neoplastic_ids, 'cell_type'] = 'Neoplastic'
write.table(cellnames_types, file = 'cellnames_types.txt', sep = '\t', quote = F)
