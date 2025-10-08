# Prepare input data for 13C-scMFA
# from GBM12, GBM38, HF2303 scRNA-seq

# Import Libraries ------------------------------------------------------------
library(Seurat)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(RColorBrewer)

# Read scRNA-seq data ---------------------------------------------------------
# Load preprocessed data saved in single_cell_analysis/GBM12_GBM38_HF2303.R
combined = readRDS('combined.rds')

# Serine model input ----------------------------------------------------------
## Define serine related genes ------------------------------------------------
list_serine_genes = c('PHGDH', 'PSAT1', 'PSPH', # glucose-derived serine synthesis
                      'SHMT1', 'SHMT2', # serine -> glycine
                      'SDS', 'SDSL', # serine -> pyruvate
                      'SRR', # L-serine -> D-serine
                      'CBS', # serine -> cysteine
                      'SPTLC1', 'SPTLC2', 'SPTLC3', # serine -> 3-dehydro sphinganine
                      'SARS1', 'SARS2', # serine -> protein serine
                      'PTDSS1', 'PTDSS2', # serine -> phosphotidyl serine
                      'SLC1A4', 'SLC1A5', 'SLC3A2', 
                      'SLC6A9', 'SLC6A14', 'SLC7A8', 'SLC7A10', 'SLC12A4', 
                      'SLC25A1', 'SLC25A15', 'SLC25A28',
                      'SLC36A1',
                      'SLC38A1', 'SLC38A2', 'SLC38A4', 'SLC38A5', 'SLC38A7') 

# adding unlabeled serine gene set
autophagy_gene_set = c('ULK1', 'ULK2', 'ATG13', 'ATG101', 'RB1CC1')

serine_model_genes = c(list_serine_genes, autophagy_gene_set)
serine_genes_human = paste('GRCh38-', serine_model_genes, sep = '')
serine_model_genes[13] = 'SARS'
serine_genes_mouse = paste('GRCm39-', 
                           stringr::str_to_sentence(serine_model_genes), 
                           sep = '')
# visualizing serine gene expression 
serine_combined = combined[c(serine_genes_human, serine_genes_mouse)]
serine_combined = SetIdent(serine_combined, value = 'human_vs_mouse')
serine_human = subset(serine_combined, idents = 'human')
serine_mouse = subset(serine_combined, idents = 'mouse')
serine_human = serine_human[serine_genes_human]
serine_mouse = serine_mouse[serine_genes_mouse]

# Rename Genes Seurat 
# map the human and mouse gene names to a common system
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], 
                              newnames = HGNC.updated[[i]]$Suggested.Symbol) { 
  # Replace gene names in different slots of a Seurat object. 
  # Run this before integration. 
  # It only changes obj@assays$RNA@counts, @data
  # DotPlot uses data slot
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
serine_human = RenameGenesSeurat(obj = serine_human, 
                                 newnames = serine_model_genes)
serine_mouse = RenameGenesSeurat(obj = serine_mouse, 
                                 newnames = serine_model_genes)
serine_hm = merge(x = serine_human, y = serine_mouse, 
                  project = 'serine_combined_gene')

# define neoplastic cells invading through brain
serine_hm[['coarse_cell_anno_3']] = as.character(serine_hm$cell_anno)
serine_hm@meta.data[serine_hm$brain_vs_tumor == 'brain' &
                      serine_hm$coarse_cell_anno_2 == 'Neoplastic', 
                    'coarse_cell_anno_3'] = 'Invading_Neoplastic'
levels(factor(serine_hm$coarse_cell_anno_3))


serine_hm$coarse_cell_anno_3 = factor(serine_hm$coarse_cell_anno_3,
                                      levels = c('Neoplastic', 
                                                 'Invading_Neoplastic',
                                                 'Neuron',
                                                 'Astrocyte', 
                                                 'Oligodendrocyte', 
                                                 'OPC', 
                                                 'Myeloid', 
                                                 'Ependymal',
                                                 'Vascular',
                                                 'Not Determined'
                                                 ))
table(serine_hm$coarse_cell_anno_3)
# Neoplastic Invading_Neoplastic            Neuron_1            Neuron_2 
# 12648                 181                2991                 279 
# Astrocyte     Oligodendrocyte                 OPC             Myeloid 
# 869                 254                  78                2043 
# Ependymal            Vascular 
# 62                 494

# Neoplastic Invading_Neoplastic              Neuron           Astrocyte 
# 12648                 181                3270                 869 
# Oligodendrocyte                 OPC             Myeloid           Ependymal 
# 254                  78                2055                  14 
# Vascular      Not Determined 
# 494                  36 

saveRDS(serine_hm, 'serine_hm.rds')

serine_hm = SetIdent(serine_hm, value = 'coarse_cell_anno_3')

p = DotPlot_scCustom(serine_hm, 
                     features = serine_model_genes,
                     idents = c('Neoplastic', 'Invading_Neoplastic',
                                'Neuron', 'Astrocyte'),
                     flip_axes = T, 
                     x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')))#, 
                         #limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_combined_serine_mouse_models.pdf', plot = p, 
       device = 'pdf', height = 10, width = 5, units = 'in')

## Create block diagonal matrix for each patient -------------------------------
# normalized gene expression of each cell type are combined in 
# a block diagonal matrix: 
# O - -
# - O -
# - - O
# Where O is the normalized gene expression matrix of a cell type with 
# the shape of number of genes by number of cells of that cell type
serine_hm = SetIdent(serine_hm, value = c('brain_vs_tumor'))
brain = subset(serine_hm, idents = 'brain')
tumor = subset(serine_hm, idents = 'tumor')

brain = SetIdent(brain, value = c('coarse_cell_anno_2'))
tumor = SetIdent(tumor, value = c('coarse_cell_anno_2'))

# Neoplastic cells in tumor samples
tumor = subset(tumor, idents = 'Neoplastic')
tumor = SetIdent(tumor, value = 'mouse_model')
table(tumor@active.ident)
# HF2303  GBM38  GBM12 
# 3706    5950   2992 

# Neoplastic cells invading through brain
tumor_invasion = subset(brain, idents = 'Neoplastic')
tumor_invasion = SetIdent(tumor_invasion, value = 'mouse_model')
table(tumor_invasion@active.ident)
# HF2303  GBM38  GBM12 
# 64      12     105 

# Astrocytes in brain
astro = subset(brain, idents = 'Astrocyte')
astro = SetIdent(astro, value = 'mouse_model')
table(astro@active.ident)
# HF2303  GBM38  GBM12 
# 180     253    406 

# Neurons in brain
neuron = subset(brain, idents = c('Neuron_1', 'Neuron_2'))
neuron = SetIdent(neuron, value = 'mouse_model')
table(neuron@active.ident)
# HF2303  GBM38  GBM12 
# 447     712    2099 

### save block diagonal matrix for each mouse model separately -----------------
for (mouse_model in unique(combined$mouse_model)){ 
  
  # Neoplastic cells in tumor samples
  tumor_mm = subset(tumor, idents = mouse_model)
  tumor_ge = GetAssayData(tumor_mm, slot = 'data')
  dim(tumor_ge) 
  dimnames(tumor_ge)[[1]] = paste(dimnames(tumor_ge)[[1]], 'g', sep = '_')
  tumor_ge = as.data.frame(tumor_ge)
  
  # Neoplastic cells invading through brain
  tumor_invasion_mm = subset(tumor_invasion, idents = mouse_model)
  tumor_invasion_ge = GetAssayData(tumor_invasion_mm, slot = 'data')
  dim(tumor_invasion_ge) 
  dimnames(tumor_invasion_ge)[[1]] = paste(dimnames(tumor_invasion_ge)[[1]], 
                                           'i', sep = '_')
  tumor_invasion_ge = as.data.frame(tumor_invasion_ge)
  
  # Astrocytes in brain
  astro_mm = subset(astro, idents = mouse_model)
  astro_ge = GetAssayData(astro_mm, slot = 'data')
  dim(astro_ge) 
  dimnames(astro_ge)[[1]] = paste(dimnames(astro_ge)[[1]], 'a', sep = '_')
  astro_ge = as.data.frame(astro_ge)
  
  # Neurons in brain
  neuron_mm = subset(neuron, idents = mouse_model)
  neuron_ge = GetAssayData(neuron_mm, slot = 'data')
  dim(neuron_ge) 
  dimnames(neuron_ge)[[1]] = paste(dimnames(neuron_ge)[[1]], 'n', sep = '_')
  neuron_ge = as.data.frame(neuron_ge)
  
  combined_scMFA = bind_rows(astro_ge, neuron_ge, tumor_invasion_ge, tumor_ge)
  combined_scMFA[is.na(combined_scMFA)] = 0
  
  write.csv(combined_scMFA, file = paste(mouse_model, '4c.csv', sep = '_'), 
            quote = F)
  
  neoplastic_tumor_ids = colnames(tumor_mm)
  neoplastic_brain_ids = colnames(tumor_invasion_mm)
  astro_brain_ids = colnames(astro_mm)
  neuron_brain_ids = colnames(neuron_mm)
  
  ids = c(astro_brain_ids, neuron_brain_ids, neoplastic_brain_ids, 
          neoplastic_tumor_ids)
  
  cellnames_types = data.frame(ids, row.names = ids)
  cellnames_types[['cell_type']] = NA
  cellnames_types[cellnames_types$ids %in% astro_brain_ids, 
                  'cell_type'] = 'Astrocyte'
  cellnames_types[cellnames_types$ids %in% neuron_brain_ids, 
                  'cell_type'] = as.character(neuron_mm$coarse_cell_anno_2)
  cellnames_types[cellnames_types$ids %in% neoplastic_brain_ids, 
                  'cell_type'] = 'Invading_Neoplastic'
  cellnames_types[cellnames_types$ids %in% neoplastic_tumor_ids, 
                  'cell_type'] = 'Neoplastic'
  write.table(cellnames_types, file = paste(mouse_model, 'cellnames_types.txt', 
                                            sep = '_'), 
              sep = '\t', quote = F)
  
}

# Purine model input ----------------------------------------------------------
## TRP tumor -------------------------------------------------------------------
### Compare purine-related genes between cell types in TRP tumor ----
mtrx = readRDS('SP3_res01_2507.rds')
list_purine_genes = c('Ppat', 'Gart', 'Pfas', 'Paics', 'Adsl', 'Atic', 
                      'Prps1l1', 'Prps1', 'Prps2', # de novo IMP
                      'Hprt', # HPX_IMP
                      'Nt5c2', 'Nt5c', 'Nt5e', 'Nt5m', 
                      'Nt5c1a', 'Nt5c1b', # IMP_out
                      'Adss', 'Adssl1', # IMP_AMP 
                      'Ampd1', 'Ampd2', 'Ampd3', # AMP_IMP
                      'Adk', # ADE_AMP
                      'Ada', 'Lacc1', # ADE_INO
                      'Pnp', # nucleoside to nucleobase
                      'Impdh1', 'Impdh2', 'Gmps', # IMP_GMP
                      'Guk1', # GMP_GDP
                      'Rrm2b', 'Rrm1', 'Rrm2', 'Nme6', 'Ak9', 'Nme7', 'Nme1', 
                      'Nme2', 'Nme3', 'Nme4', # GDP_out
                      'Ak1', 'Ak2', 'Ak3', 'Ak4', 'Ak5', 'Ak6', 'Ak7', 'Ak8', # AMP_out
                      'Slc28a1', 'Slc28a2', 'Slc28a3', # CNT
                      'Slc29a1', 'Slc29a2', 'Slc29a3', # ENT
                      'Slc43a3') # ENBT1

mtrx = SetIdent(mtrx, value = 'cell_anno_2')
mtrx$cell_anno_2 = factor(as.character(mtrx$cell_anno_2), 
                         levels = c('Neoplastic',
                                    'Myeloid', 
                                    "Oligodendrocyte",
                                    "OPC",
                                    "Lymphoid",
                                    "Vascular"))

p = DotPlot_scCustom(mtrx, features = list_purine_genes, 
                     flip_axes = T, 
                     x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_TRP_SP3_purine_input_scFBA_MFA_2507.pdf', plot = p, 
       device = 'pdf', height = 12, width = 4.5, units = 'in')


### Create block diagonal matrix for 13C-scMFA --------------------------------
# normalized gene expression of each cell type are combined in 
# a block diagonal matrix: 
# O - -
# - O -
# - - O
# Where O is the normalized gene expression matrix of a cell type with 
# the shape of number of genes by number of cells of that cell type

list_purine_human_genes = c('PPAT', 'GART', 'PFAS', 'PAICS', 'ADSL', 'ATIC', 
                      'PRPS1L1', 'PRPS1', 'PRPS2', # de novo IMP
                      'HPRT1', # HPX_IMP
                      'NT5C2', 'NT5C', 'NT5E', 'NT5M', 
                      'NT5C1A', 'NT5C1B', # IMP_out
                      'ADSS', 'ADSSL1', # IMP_AMP 
                      'AMPD1', 'AMPD2', 'AMPD3', # AMP_IMP
                      'ADK', # ADE_AMP
                      'ADA', 'LACC1', # ADE_INO
                      'PNP', # nucleoside to nucleobase
                      'IMPDH1', 'IMPDH2', 'GMPS', # IMP_GMP
                      'GUK1', # GMP_GDP
                      'RRM2B', 'RRM1', 'RRM2', 'NME6', 'AK9', 'NME7', 'NME1', 
                      'NME2', 'NME3', 'NME4', # GDP_out
                      'AK1', 'AK2', 'AK3', 'AK4', 'AK5', 'AK6', 'AK7', 'AK8', # AMP_out
                      'SLC28A1', 'SLC28A2', 'SLC28A3', # CNT
                      'SLC29A1', 'SLC29A2', 'SLC29A3', # ENT
                      'SLC43A3') # ENBT1

mtrx = SetIdent(mtrx, value = 'cell_anno_2')

# malignant cells
mal = subset(mtrx, idents = 'Neoplastic')

mal_ge = GetAssayData(mal, slot = 'data') # normalized
mal_ge = mal_ge[list_purine_genes, ]
dim(mal_ge) #54 2670
dimnames(mal_ge)[[1]] = paste(list_purine_human_genes, 'g', sep = '_')
mal_ge = as.data.frame(mal_ge)

# myeloid cells
myeloid = subset(mtrx, idents = 'Myeloid')

my_ge = GetAssayData(myeloid, slot = 'data') # normalized
my_ge = my_ge[list_purine_genes, ]
dim(my_ge) #54 3971
dimnames(my_ge)[[1]] = paste(list_purine_human_genes, 'm', sep = '_')
my_ge = as.data.frame(my_ge)

combined = bind_rows(my_ge, mal_ge)
combined[is.na(combined)] = 0

myeloid_ids = colnames(my_ge)
neoplastic_ids = colnames(mal_ge)
ids = c(myeloid_ids, neoplastic_ids)

# save cell IDs and cell types
cellnames_types = data.frame(ids, row.names = ids)
cellnames_types[['cell_type']] = NA

cellnames_types[cellnames_types$ids %in% myeloid_ids, 'cell_type'] = 'Myeloid'
cellnames_types[cellnames_types$ids %in% neoplastic_ids, 'cell_type'] = 'Neoplastic'
write.table(cellnames_types, file = paste('SP3', 'purine_cellnames_types.txt', 
                                          sep = '_'), sep = '\t', quote = F)
# save block diagonal matrix of gene expression  
write.csv(combined, file = paste("SP3", 'purine_scRNA.csv', sep = '_'), 
          quote = F)


## GBM38 tumor ----------------------------------------------------------------
setwd("./Baharan Meghdadi/Baharan-Deepak/ML_manuscript/codes/single_cell_analysis/GBM12_GBM38_HF2303_figs/human_neoplastic_mouse_non_neoplastic/updated_cl210/cl13/")
three_mice = readRDS('combined_res01.rds')
three_mice = SetIdent(three_mice, value = 'orig.ident')

mtrx = subset(three_mice, idents = 'GBM38_tumor_11913_SP1')

list_purine_human_genes_2 = c('PPAT', 'GART', 'PFAS', 'PAICS', 'ADSL', 'ATIC', 
                            'PRPS1L1', 'PRPS1', 'PRPS2', # de novo IMP
                            'HPRT1', # HPX_IMP
                            'NT5C2', 'NT5C', 'NT5E', 'NT5M', 
                            'NT5C1A', 'NT5C1B', # IMP_out
                            'ADSS2', 'ADSS1', # IMP_AMP 
                            'AMPD1', 'AMPD2', 'AMPD3', # AMP_IMP
                            'ADK', # ADE_AMP
                            'ADA', 'LACC1', # ADE_INO
                            'PNP', # nucleoside to nucleobase
                            'IMPDH1', 'IMPDH2', 'GMPS', # IMP_GMP
                            'GUK1', # GMP_GDP
                            'RRM2B', 'RRM1', 'RRM2', 'NME6', 'AK9', 'NME7', 'NME1', 
                            'NME2', 'NME3', 'NME4', # GDP_out
                            'AK1', 'AK2', 'AK3', 'AK4', 'AK5', 'AK6', 'AK7', 'AK8', # AMP_out
                            'SLC28A1', 'SLC28A2', 'SLC28A3', # CNT
                            'SLC29A1', 'SLC29A2', 'SLC29A3', # ENT
                            'SLC43A3') # ENBT1

purine_genes_human = paste('GRCh38-', list_purine_human_genes_2, sep = '')
purine_genes_mouse = paste('GRCm39-', 
                           stringr::str_to_sentence(list_purine_human_genes_2), 
                           sep = '')
purine_genes_mouse[10] = "GRCm39-Hprt"
purine_genes_mouse[17] = "GRCm39-Adss"
purine_genes_mouse[18] = "GRCm39-Adssl1"
# visualizing purine gene expression 
purine_combined = mtrx[c(purine_genes_human, purine_genes_mouse)]
purine_combined = SetIdent(purine_combined, value = 'human_vs_mouse')
purine_human = subset(purine_combined, idents = 'human')
purine_mouse = subset(purine_combined, idents = 'mouse')
purine_human = purine_human[purine_genes_human]
purine_mouse = purine_mouse[purine_genes_mouse]

list_purine_human_genes_3 = c('PPAT', 'GART', 'PFAS', 'PAICS', 'ADSL', 'ATIC', 
                              'PRPS1L1', 'PRPS1', 'PRPS2', # de novo IMP
                              'HPRT1', # HPX_IMP
                              'NT5C2', 'NT5C', 'NT5E', 'NT5M', 
                              'NT5C1A', 'NT5C1B', # IMP_out
                              'ADSS', 'ADSSL1', # IMP_AMP 
                              'AMPD1', 'AMPD2', 'AMPD3', # AMP_IMP
                              'ADK', # ADE_AMP
                              'ADA', 'LACC1', # ADE_INO
                              'PNP', # nucleoside to nucleobase
                              'IMPDH1', 'IMPDH2', 'GMPS', # IMP_GMP
                              'GUK1', # GMP_GDP
                              'RRM2B', 'RRM1', 'RRM2', 'NME6', 'AK9', 'NME7', 'NME1', 
                              'NME2', 'NME3', 'NME4', # GDP_out
                              'AK1', 'AK2', 'AK3', 'AK4', 'AK5', 'AK6', 'AK7', 'AK8', # AMP_out
                              'SLC28A1', 'SLC28A2', 'SLC28A3', # CNT
                              'SLC29A1', 'SLC29A2', 'SLC29A3', # ENT
                              'SLC43A3') # ENBT1
# Rename Genes Seurat 
# map the human and mouse gene names to a common system
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], 
                              newnames = HGNC.updated[[i]]$Suggested.Symbol) { 
  # Replace gene names in different slots of a Seurat object. 
  # Run this before integration. 
  # It only changes obj@assays$RNA@counts, @data
  # DotPlot uses data slot
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
purine_human = RenameGenesSeurat(obj = purine_human, 
                                 newnames = list_purine_human_genes_3)
purine_mouse = RenameGenesSeurat(obj = purine_mouse, 
                                 newnames = list_purine_human_genes_3)
purine_hm = merge(x = purine_human, y = purine_mouse, 
                  project = 'GBM38_SP1_purine_combined_gene')


saveRDS(purine_hm, 'purine_hm.rds')

purine_hm = SetIdent(purine_hm, value = 'cell_anno')
purine_hm$cell_anno = factor(as.character(purine_hm$cell_anno), 
                             levels = c("Neoplastic",
                                        "Myeloid",
                                        "Oligodendrocyte",
                                        "OPC",
                                        "Astrocyte",
                                        "Neuron",
                                        "Vascular"))
### Compare purine-related genes between cell types ----
p = DotPlot_scCustom(purine_hm, 
                     features = list_purine_human_genes_3,
                     flip_axes = T, 
                     x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')),
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_GBM38_SP1_purine_250729.pdf', plot = p, 
       device = 'pdf', height = 12, width = 4.5, units = 'in')

### Create block diagonal matrix for 13C-scMFA --------------------------------
# malignant cells
mal = subset(purine_hm, idents = 'Neoplastic')

mal_ge = GetAssayData(mal, slot = 'data') # normalized

dim(mal_ge) #54 5950
dimnames(mal_ge)[[1]] = paste(list_purine_human_genes_3, 'g', sep = '_')
mal_ge = as.data.frame(mal_ge)

# myeloid cells
myeloid = subset(purine_hm, idents = 'Myeloid')

my_ge = GetAssayData(myeloid, slot = 'data') # normalized

dim(my_ge) #54 858
dimnames(my_ge)[[1]] = paste(list_purine_human_genes_3, 'm', sep = '_')
my_ge = as.data.frame(my_ge)

combined = bind_rows(my_ge, mal_ge)
combined[is.na(combined)] = 0

myeloid_ids = colnames(my_ge)
neoplastic_ids = colnames(mal_ge)
ids = c(myeloid_ids, neoplastic_ids)

# save cell IDs and cell types
cellnames_types = data.frame(ids, row.names = ids)
cellnames_types[['cell_type']] = NA

cellnames_types[cellnames_types$ids %in% myeloid_ids, 'cell_type'] = 'Myeloid'
cellnames_types[cellnames_types$ids %in% neoplastic_ids, 'cell_type'] = 'Neoplastic'
write.table(cellnames_types, file = paste('SP1', 'purine_cellnames_types.txt', 
                                          sep = '_'), sep = '\t', quote = F)
# save block diagonal matrix of gene expression for each patient
write.csv(combined, file = paste("SP1", 'purine_scRNA.csv', sep = '_'), 
          quote = F)
