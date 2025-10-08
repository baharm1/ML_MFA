# Calculate Metabolite Accumulation Score (MAS) for serine
# scRNA-seq dataset from Darmanis et al., Cell Reports 2017

# Import Libraries ------------------------------------------------------------
library(Seurat) # version 4.2.0
library(scCustomize) # version 1.0.0
library(ggplot2) # version 3.4.2
library(Rmagic) # version 2.0.3
library(RColorBrewer) # version 1.1-3

# Read Darmanis dataset -------------------------------------------------------
# Load preprocessed data saved in scRNA_analysis_Darmanis.R
seuobj = readRDS('darmanis.rds')

# Define serine related genes -------------------------------------------------
# Serine synthesis pathway
serine_prod = c('PHGDH', 'PSAT1', 'PSPH') # PG to serine
# Serine catabolic enzymes
serine_cons = c('SHMT1', 'SHMT2', # serine to glycine
                'SDS', 'SDSL', # serine to pyruvate, serine to 2-aminoarcrylate
                'CBS', # serine and L-homocysteine to cystathionine
                'SPTLC1', 'SPTLC2', 'SPTLC3', # serine and palmitoyl coa to 3-dehydro sphinganine
                'SARS', 'SARS2', # seryl t-RNA synthetase
                'SRR', # L-serine to D-serine
                'PTDSS1', 'PTDSS2') # serine to phosphatidyl serine
# Serine transporters
serine_tran = c('SLC1A4', 'SLC1A5', 'SLC3A2', 
                'SLC6A9', 'SLC6A14', 'SLC7A8', 'SLC7A10', 'SLC12A4', 
                'SLC25A1', 'SLC25A15', 'SLC25A28',
                'SLC36A1',
                'SLC38A1', 'SLC38A2', 'SLC38A4', 'SLC38A5', 'SLC38A7')

# Calculate metabolite accumulation score (MAS) -------------------------------
# Impute normalized gene expression
magic_obj = magic(t(GetAssayData(seuobj, slot = 'data')))
magic_res = magic_obj$result

# Create scores for serine synthesis, consumption, and exchange
serine_prod_magic = rowSums(magic_res[, serine_prod])
serine_cons_magic = rowSums(magic_res[, serine_cons])
serine_tran_magic = rowSums(magic_res[, serine_tran])

# MAS 
serine_prod_cons = serine_prod_magic - serine_cons_magic

# Add scores to metadata
seuobj@meta.data[['serine_prod_cons']] = serine_prod_cons
seuobj@meta.data[['serine_prod']] = serine_prod_magic
seuobj@meta.data[['serine_cons']] = serine_cons_magic
seuobj@meta.data[['serine_tran']] = serine_tran_magic

# Show distribution of scores for various cell types
p = VlnPlot(seuobj, features = c('serine_tran'), 
            group.by = 'major_cell_type') + 
  NoLegend() +
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white")

ggsave(filename = 'vlnplot_major_cell_type_serine_tran.pdf', plot = p, 
       device = 'pdf', height = 6, width = 6, units = 'in')

p = VlnPlot(seuobj, features = c('serine_prod_cons'), 
            group.by = 'Location') + 
  NoLegend() +
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white") 

ggsave(filename = 'vlnplot_location_serine_prod_cons.pdf', plot = p, 
       device = 'pdf', height = 6, width = 3, units = 'in')

# Show scaled gene expression of serine related genes for various cell types
p = DotPlot_scCustom(seuobj, features = c(serine_prod, serine_cons, serine_tran), 
                     group.by = 'major_cell_type',
                     flip_axes = T, x_lab_rotate = T) +
  theme(text = element_text(size = 14)) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))

ggsave(filename = 'dotplot_major_cell_type_serine_all.pdf', plot = p, 
       device = 'pdf', height = 10, width = 6, units = 'in')

saveRDS(seuobj, 'darmanis.rds')
