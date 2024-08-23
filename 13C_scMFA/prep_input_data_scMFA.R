# Prepare input data for 13C-scMFA

# Import Libraries ------------------------------------------------------------
library(Seurat)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(RColorBrewer)

# Read scRNA-seq data ---------------------------------------------------------
# Load preprocessed data saved in scRNA_analysis_Darmanis.R
# Load preprocessed data saved in scRNA_analysis_our_data.R
seuobj = readRDS('darmanis.rds')
seuobj = SetIdent(seuobj, value = c('major_cell_type'))
gbm = readRDS('gbm_w_metadata.rds')
gbm = SetIdent(gbm, value = 'cell_anno_2')

# Serine model input ----------------------------------------------------------
## Define serine related genes ------------------------------------------------
list_serine_genes = c('PHGDH', 'PSAT1', 'PSPH', # glucose-derived serine synthesis
                      'SHMT1', 'SHMT2', # serine -> glycine
                      'SDS', 'SDSL', # serine -> pyruvate
                      'SRR', # L-serine -> D-serine
                      'CBS', # serine -> cysteine
                      'SPTLC1', 'SPTLC2', 'SPTLC3', # serine -> 3-dehydro sphinganine
                      'SARS', 'SARS2', # serine -> protein serine
                      'PTDSS1', 'PTDSS2', # serine -> phosphotidyl serine
                      'SLC1A4', 'SLC1A5', 'SLC3A2', 
                      'SLC6A9', 'SLC6A14', 'SLC7A8', 'SLC7A10', 'SLC12A4', 
                      'SLC25A1', 'SLC25A15', 'SLC25A28',
                      'SLC36A1',
                      'SLC38A1', 'SLC38A2', 'SLC38A4', 'SLC38A5', 'SLC38A7')

# adding unlabeled serine gene set
autophagy_gene_set = c('ULK1', 'ULK2', 'ATG13', 'C12orf44', 'RB1CC1')
autophagy_gene_set3 = c('ULK1', 'ULK2', 'ATG13', 'ATG101', 'RB1CC1')

list_serine_genes2 = c(list_serine_genes, autophagy_gene_set)
list_serine_genes3 = c(list_serine_genes, autophagy_gene_set3)

## Create block diagonal matrix for each patient -------------------------------
# normalized gene expression of each cell type are combined in 
# a block diagonal matrix: 
# O - -
# - O -
# - - O
# Where O is the normalized gene expression matrix of a cell type with 
# the shape of number of genes by number of cells of that cell type

# Neoplastic cells are from our patient cohort
gbm@meta.data[['patientSite']] = paste(gbm$patient, gbm$site, sep = '_')
table(gbm$patientSite)
mal = subset(gbm, idents = c('AC-like', 'NPC-like', 'MES-like', 'OPC-like', 'CellCycle'))
mal = SetIdent(mal, value = 'patientSite')

# Astrocytes and neurons are from Darmanis dataset
astro = subset(seuobj, idents = c('Astrocyte'))
astro_ge = GetAssayData(astro, slot = 'data')
astro_ge = astro_ge[list_serine_genes2, ]
dim(astro_ge) # 36 * 88

dimnames(astro_ge)[[1]] = paste(dimnames(astro_ge)[[1]], 'a', sep = '_')
astro_ge = as.data.frame(astro_ge)

neuron = subset(seuobj, idents = c('Neuron'))
neuron_ge = GetAssayData(neuron, slot = 'data')
neuron_ge = neuron_ge[list_serine_genes2, ]
dim(neuron_ge) # 36 * 21

dimnames(neuron_ge)[[1]] = paste(dimnames(neuron_ge)[[1]], 'n', sep = '_')
neuron_ge = as.data.frame(neuron_ge)

metadata = seuobj@meta.data

astrocyte_ids = rownames(metadata[metadata$major_cell_type == 'Astrocyte', ])
neuron_ids = rownames(metadata[metadata$major_cell_type == 'Neuron', ])

# save block diagonal matrix for each patient separately 
for (patient_site in unique(mal$patientSite)){ 
  
  mal_ps = subset(mal, idents = patient_site)
  mal_ge = GetAssayData(mal_ps, slot = 'data')
  mal_ge = mal_ge[list_serine_genes3, ]
  dim(mal_ge) 
  
  dimnames(mal_ge)[[1]] = paste(dimnames(mal_ge)[[1]], 'g', sep = '_')
  mal_ge = as.data.frame(mal_ge)
  
  combined = bind_rows(astro_ge, neuron_ge, mal_ge)
  combined[is.na(combined)] = 0
  
  write.csv(combined, file = paste(patient_site, '3c_scRNA.csv', sep = '_'), 
            quote = F)
  
  neoplastic_ids = colnames(mal_ps)
  
  ids = c(astrocyte_ids, neuron_ids, neoplastic_ids)
  
  cellnames_types = data.frame(ids, row.names = ids)
  cellnames_types[['cell_type']] = NA
  cellnames_types[cellnames_types$ids %in% astrocyte_ids, 'cell_type'] = 'Astrocyte'
  cellnames_types[cellnames_types$ids %in% neuron_ids, 'cell_type'] = 'Neuron'
  cellnames_types[cellnames_types$ids %in% neoplastic_ids, 'cell_type'] = 'Neoplastic'
  write.table(cellnames_types, file = paste(patient_site, 'cellnames_types.txt', 
                                            sep = '_'), sep = '\t', quote = F)
  
}

# Purine model input ----------------------------------------------------------
## Define purine related genes ------------------------------------------------
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

## Identify cell types that should be included in the 13C-scMFA model ----
df = data.frame(sample = gbm@assays$RNA@data@Dimnames[[2]], 
                anno = as.character(gbm$cell_anno_2))
table(gbm$anno)
df[df$anno == 'AC-like', 'anno'] = 'Neoplastic'
df[df$anno == 'NPC-like', 'anno'] = 'Neoplastic'
df[df$anno == 'OPC-like', 'anno'] = 'Neoplastic'
df[df$anno == 'MES-like', 'anno'] = 'Neoplastic'
df[df$anno == 'CellCycle', 'anno'] = 'Neoplastic'
rownames(df) = df$sample

gbm = AddMetaData(gbm, metadata = df)                      

### Compare purine-related genes between cell types ----
# Extended Fig. 15A
p = DotPlot_scCustom(gbm, features = list_purine_genes, 
                     group.by = 'anno',
                     flip_axes = T, 
                     x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_scwna_purine_input_scFBA_MFA_2407.pdf', plot = p, 
       device = 'pdf', height = 10, width = 4.5, units = 'in')

### Create module scores for purine reactions ----
# Extended Fig. 15B
gbm = AddModuleScore(gbm, features = list(c('PPAT', 'GART', 
                                            'PFAS', 'PAICS', 
                                            'ADSL', 'ATIC', 
                                            'PRPS1L1', 'PRPS1', 
                                            'PRPS2')), name = 'denovo_IMP')

gbm = AddModuleScore(gbm, features = list(c('NT5C2', 'NT5C', 
                                            'NT5E', 'NT5DC4', 'NT5M', 
                                            'NT5C1A', 'NT5C1B')), 
                     name = 'nucleotide_nucleoside')

gbm = AddModuleScore(gbm, features = list(c('IMPDH1', 'IMPDH2', 'GMPS')), 
                     name = 'IMP_GMP')

gbm = AddModuleScore(gbm, features = list(c('AMPD1', 'AMPD2', 'AMPD3')), 
                     name = 'AMP_IMP')

gbm = AddModuleScore(gbm, features = list(c('ADSS', 'ADSSL1', 'ADSL')), 
                     name = 'IMP_AMP')

gbm = AddModuleScore(gbm, features = list(c('ADA', 'ADA2', 'LACC1')), 
                     name = 'ADE_INO')

gbm = AddModuleScore(gbm, features = list(c('PNP', 'LACC1')), 
                     name = 'nucleoside_nucleobase')

gbm = AddModuleScore(gbm, features = list(c('AK1', 'AK2', 'AK3', 
                                            'AK4', 'AK5', 'AK6', 
                                            'AK7', 'AK8', 'AK9')), 
                     name = 'AMP_ADP')

gbm = AddModuleScore(gbm, features = list(c('SLC28A1', 'SLC28A2', 'SLC28A3', # CNT
                                            'SLC29A1', 'SLC29A2', 'SLC29A3', # ENT
                                            'SLC43A3')), 
                     name = 'nucleobase_transporter')

p = DotPlot_scCustom(gbm, features = c('HPRT1', 
                                       'denovo_IMP1', 
                                       'IMP_GMP1',
                                       'nucleotide_nucleoside1',
                                       'nucleoside_nucleobase1',
                                       'nucleobase_transporter1'), 
                     group.by = 'anno',
                     x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_scwna_purine_modulescore_2407.pdf', plot = p, 
       device = 'pdf', height = 4, width = 6, units = 'in')


# Compare purine-related genes between cell types patient-wise
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

## Create block diagonal matrix for each patient ------------------------------
# normalized gene expression of each cell type are combined in 
# a block diagonal matrix: 
# O - -
# - O -
# - - O
# Where O is the normalized gene expression matrix of a cell type with 
# the shape of number of genes by number of cells of that cell type

# patinet tumors: enhancing and nonenhancing
gbm@meta.data[['patientSite']] = paste(gbm$patient, gbm$site, sep = '_')
table(gbm$patientSite)
# malignant cells
mal = subset(gbm, idents = c('AC-like', 'NPC-like', 'MES-like', 'OPC-like', 'CellCycle'))
mal = SetIdent(mal, value = 'patientSite')
# myeloid cells
myeloid = subset(gbm, idents = c('Myeloid'))
myeloid = SetIdent(myeloid, value = 'patientSite')

for (patient_site in unique(mal$patientSite)){ 
  
  # malignant cells
  mal_ps = subset(mal, idents = patient_site)
  mal_ge = GetAssayData(mal_ps, slot = 'data') # normalized
  mal_ge = mal_ge[list_purine_genes, ]
  dim(mal_ge) 
  dimnames(mal_ge)[[1]] = paste(dimnames(mal_ge)[[1]], 'g', sep = '_')
  mal_ge = as.data.frame(mal_ge)
  
  # myeloid cells
  my_ps = subset(myeloid, idents = patient_site)
  my_ge = GetAssayData(my_ps, slot = 'data') # normalized
  my_ge = my_ge[list_purine_genes, ]
  dim(my_ge) 
  dimnames(my_ge)[[1]] = paste(dimnames(my_ge)[[1]], 'm', sep = '_')
  my_ge = as.data.frame(my_ge)
  
  combined = bind_rows(my_ge, mal_ge)
  combined[is.na(combined)] = 0
  
  myeloid_ids = colnames(my_ps)
  neoplastic_ids = colnames(mal_ps)
  ids = c(myeloid_ids, neoplastic_ids)
  
  # save cell IDs and cell types
  cellnames_types = data.frame(ids, row.names = ids)
  cellnames_types[['cell_type']] = NA
  
  cellnames_types[cellnames_types$ids %in% myeloid_ids, 'cell_type'] = 'Myeloid'
  cellnames_types[cellnames_types$ids %in% neoplastic_ids, 'cell_type'] = 'Neoplastic'
  write.table(cellnames_types, file = paste(patient_site, 'purine_cellnames_types.txt', 
                                            sep = '_'), sep = '\t', quote = F)
  # save block diagonal matrix of gene expression for each patient
  write.csv(combined, file = paste(patient_site, 'purine_scRNA.csv', sep = '_'), 
            quote = F)
  
}
