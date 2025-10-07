# Preprocess scRNA-seq dataset from Darmanis et al., Cell Reports 2017

# Import Libraries ------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(reshape2)
library(Polychrome)  
library(scales)
library(DropletUtils)

# Plot Functions --------------------------------------------------------------
# Proportion plot
plot_integrated_clusters = function (srat, cluster_x, cluster_y, is_cluster_x_int, 
                                     filename, height = 5.5, width = 12) { 
  # take an integrated Seurat object, plot distributions over orig.ident
  
  count_table <- table(srat@meta.data[[cluster_x]], srat@meta.data[[cluster_y]])
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)
  data("palette36")
  
  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  if (is_cluster_x_int){
    sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)), decreasing = T))
    cluster_size$cluster <- factor(cluster_size$cluster, levels = sorted_labels)
  } else {
    cluster_size$cluster = factor(cluster_size$cluster)
    sorted_labels = levels(cluster_size$cluster)
  }
  
  melt_mtx$cluster <- factor(melt_mtx$cluster, levels = sorted_labels)
  
  colnames(melt_mtx)[2] <- "dataset"
  
  
  p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + 
    geom_bar(position="dodge", stat="identity", fill = "grey60") + 
    theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("")
  p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
    geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
    ylab("Fraction of cells in each dataset") + xlab(cluster_x) + 
    theme(legend.position="top")
  
  p = p2 + p1 + plot_layout(widths = c(3, 1))
  ggsave(filename = paste(filename, '.pdf', sep = ''),
         plot = p, device = 'pdf', width = width, height = height, units = 'in')
  
}

# Read scRNA-seq data ---------------------------------------------------------
# Read raw counts
read_counts = read.delim(file = 'GBM_data_and_metadata/GBM_raw_gene_counts.csv',
                         header = T, row.names = 1, sep = ' ') # 23368 * 3589
# Read metadata
metadata = read.delim(file = 'GBM_data_and_metadata/GBM_metadata.csv',
                      header = T, row.names = 1, sep = ' ')
rownames(metadata) = paste('X', rownames(metadata), sep = '')

table(metadata$Selection)
table(metadata$Location)
table(metadata$Sample.name)

seuobj = CreateSeuratObject(read_counts, project = "Darmanis2017", 
                            meta.data = metadata)

# Preprocessing steps based on Seurat workflow --------------------------------
seuobj = NormalizeData(seuobj)
seuobj = FindVariableFeatures(seuobj)
seuobj = ScaleData(seuobj, features = rownames(seuobj))
seuobj = RunPCA(seuobj)
ElbowPlot(seuobj)

seuobj = FindNeighbors(seuobj, dims = 1:20)
seuobj = FindClusters(seuobj)
seuobj = RunUMAP(seuobj, dims = 1:20)
DimPlot(seuobj, reduction = 'umap')

seuobj@meta.data[seuobj@meta.data$Location == 'Distant', 'Location'] = 'Periphery'

p = DimPlot(seuobj, reduction = 'umap', group.by = 'Location')
ggsave(filename = 'umap_Location.pdf', plot = p, device = 'pdf', 
       height = 6, width = 8, units = 'in')

# Cell annotation -------------------------------------------------------------
DE_markers = FindAllMarkers(seuobj, only.pos = T)

DE_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

p = DoHeatmap(seuobj, features = top10$gene) + NoLegend()
ggsave(filename = 'heatmap_DE_markers.pdf', plot = p, device = 'pdf', 
       height = 20, width = 6, units = 'in')

unique(metadata$Cluster_2d)
unique(metadata$Cluster_2d_color)
darmanis_colors = c('11' = "#8C564B", '2' = "#AEC7E8", '9' = "#9467BD", 
                    '8' = "#FF9896", '7' = "#D62728", '5' = "#2CA02C", 
                    '3' = "#FF7F0E", '1' = "#1F77B4", '10' = "#C5B0D5", 
                    '6' = "#98DF8A", '12' = "#C49C94", '4' = "#FFBB78")

darmanis_colors_2 = darmanis_colors[order(as.integer(names(darmanis_colors)))]
scales::show_col(darmanis_colors_2)

p = DimPlot(seuobj, reduction = 'umap', group.by = 'Cluster_2d', 
            cols = darmanis_colors_2, label = T, repel = T)
ggsave(filename = 'umap_Cluster_2d.pdf', plot = p, device = 'pdf', 
       height = 6, width = 6, units = 'in')

seuobj@meta.data[['cell_type']] = NA

seuobj@meta.data[seuobj@meta.data$Cluster_2d == '1', 'cell_type'] = 'Neoplastic cells 1'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '2', 'cell_type'] = 'Oligodendrocytes'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '3', 'cell_type'] = 'Vascular cells 1'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '4', 'cell_type'] = 'Neoplastic cells 2'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '5', 'cell_type'] = 'Neurons'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '6', 'cell_type'] = 'Vascular cells 2'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '7', 'cell_type'] = 'Myeloid cells 1'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '8', 'cell_type'] = 'Myeloid cells 2'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '9', 'cell_type'] = 'OPCs'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '10', 'cell_type'] = 'Astrocytes'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '11', 'cell_type'] = 'Neoplastic cells 3'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '12', 'cell_type'] = 'Vascular cells 3'

darmanis_colors = c('Neoplastic cells 1' = "#1F77B4", 'Neoplastic cells 2' = "#FFBB78",
                    'Neoplastic cells 3' = "#8C564B", 'Vascular cells 1' = "#FF7F0E",
                    'Vascular cells 2' = "#98DF8A", 'Vascular cells 3' = "#C49C94",
                    'Myeloid cells 1' = "#D62728", 'Myeloid cells 2' = "#FF9896",
                    'Neurons' = "#2CA02C", 'Oligodendrocytes' = "#AEC7E8",
                    'OPCs' = "#9467BD", 'Astrocytes' = "#C5B0D5")

p = DimPlot(seuobj, reduction = 'umap', group.by = 'cell_type', 
            cols = darmanis_colors, label = T, repel = T)
ggsave(filename = 'umap_cell_type.pdf', plot = p, device = 'pdf', 
       height = 6, width = 8, units = 'in')

seuobj@meta.data[['major_cell_type']] = NA
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '1', 'major_cell_type'] = 'Neoplastic'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '2', 'major_cell_type'] = 'Oligo'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '3', 'major_cell_type'] = 'Vascular'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '4', 'major_cell_type'] = 'Neoplastic'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '5', 'major_cell_type'] = 'Neuron'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '6', 'major_cell_type'] = 'Vascular'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '7', 'major_cell_type'] = 'Myeloid'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '8', 'major_cell_type'] = 'Myeloid'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '9', 'major_cell_type'] = 'OPC'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '10', 'major_cell_type'] = 'Astrocyte'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '11', 'major_cell_type'] = 'Neoplastic'
seuobj@meta.data[seuobj@meta.data$Cluster_2d == '12', 'major_cell_type'] = 'Vascular'

p = DimPlot(seuobj, reduction = 'umap', group.by = 'major_cell_type', 
            label = T, repel = T)
ggsave(filename = 'umap_major_cell_type.pdf', plot = p, device = 'pdf', 
       height = 6, width = 8, units = 'in')

plot_integrated_clusters(seuobj, 'major_cell_type', 'Location', F, 
                         'major_cell_type_Location', height = 4, width = 8)
# save Seurat object
saveRDS(seuobj, 'darmanis.rds')

# Save input data for MEBOCOST analysis ---------------------------------------
# log-normalized with scale factor of 10,000
ge_counts = GetAssayData(seuobj, slot = 'data') 
write.table(seuobj@meta.data, 
            'metadata_darmanis.csv',
            row.names = T, col.names = T, quote = F, sep = ',')

# save normalized gene expression as .h5
write10xCounts('counts_darmanis.h5',
               ge_counts, barcodes = colnames(seuobj),
               gene.id = rownames(seuobj),
               gene.symbol = rownames(seuobj),
               gene.type = "Gene Expression",
               overwrite = TRUE,
               type = c("HDF5"),
               genome = "unknown",
               version = c("3"),
               chemistry = NULL,
               original.gem.groups = NULL,
               library.ids = NULL)
