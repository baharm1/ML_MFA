# Preprocess our scRNA-seq data

# Import Libraries ------------------------------------------------------------
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(Polychrome)
library(scCustomize)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(harmony)
library(cowplot)

# Directories -----------------------------------------------------------------
directory = "./Baharan Meghdadi/Baharan-Deepak/ML_manuscript/codes/single_cell_analysis/filtered_feature_bc_matrix_GBM12_GBM38_HF2303"
fig_directory = './Baharan Meghdadi/Baharan-Deepak/ML_manuscript/codes/single_cell_analysis/GBM12_GBM38_HF2303_figs'

# Preprocessing Steps and Quality Control -------------------------------------

# PDX samples -----------------------------------------------------------------
# include both human and mouse cells
# mapped to combined human and mouse genome

setwd(directory)
sample_dir = list.files(path = '.', 
                        all.files = F, full.names = T)
sample_names = unlist(strsplit(sample_dir, '/'))[c(FALSE, TRUE)]

## read filtered count matrices
mtrx_list = c()
for (sample_no in 1:length(sample_dir)){
barcodes_dir = list.files(path = sample_dir[sample_no], 
                          pattern = 'barcodes.tsv.gz', 
                          all.files = T,
                          full.names = T)
barcodes = read.table(gzfile(barcodes_dir), stringsAsFactors = F)

matrix_dir = list.files(path = sample_dir[sample_no], 
                        pattern = 'matrix.mtx.gz', 
                        all.files = T,
                        full.names = T)
mtrx = readMM(gzfile(matrix_dir))

features_dir = list.files(path = sample_dir[sample_no], 
                          pattern = 'features.tsv.gz', 
                          all.files = T,
                          full.names = T)
features = read.table(gzfile(features_dir), stringsAsFactors = F)

dimnames(mtrx) = list(features$V2, 
                      paste(sample_names[sample_no], barcodes$V1, sep = '-'))
mtrx = CreateSeuratObject(counts = mtrx, project = sample_names[sample_no])
mtrx$orig.ident = sample_names[sample_no]
mtrx_list = append(mtrx, mtrx_list)
rm(mtrx, barcodes, features)
}

## merge filtered counts into a seurat object ----
gbm = merge(mtrx_list[[1]], y = mtrx_list[2:length(sample_names)],
            project = 'GBM12_GBM38_HF2303')

# PDX samples
gbm$orig.ident = factor(gbm$orig.ident)
gbm@active.ident = gbm$orig.ident

# add mitochondrial transcripts and mouse and human read ratios to metadata ----
gbm[['percent_mito']] = PercentageFeatureSet(gbm, 
                                             pattern = "^GRCh38-MT-|^GRCm39-mt-")

gbm[['percent_GRCm39']] = PercentageFeatureSet(gbm, pattern = "^GRCm39-")
gbm[['percent_GRCh38']] = PercentageFeatureSet(gbm, pattern = "^GRCh38-")

gbm[['GRCm39_exp']] = gbm$percent_GRCm39 * gbm$nCount_RNA / 100
gbm[['GRCh38_exp']] = gbm$percent_GRCh38 * gbm$nCount_RNA / 100

gbm = SetIdent(gbm, value = gbm@project.name)

setwd(fig_directory)
pdf(file = 'qc_vlnPlot_before.pdf', width = 6, height = 6)
VlnPlot(gbm, 
        features = c('nCount_RNA', 'nFeature_RNA', 'percent_mito'), 
        raster = F, ncol = 4, pt.size = 0)
dev.off()

# define thresholds for number of UMI counts, number of detected genes, 
# percentages of mitochondrial transcripts
pdf(file = 'qc_vlnPlot_thresholds_1.pdf', width = 10, height = 6)
p1 <- scCustomize::QC_Plots_Genes(seurat_object = gbm, low_cutoff = 500, 
                                  high_cutoff = 10000, pt.size = 0)
p2 <- QC_Plots_UMIs(seurat_object = gbm, low_cutoff = 1000, high_cutoff = 100000,
                    pt.size = 0)
p3 <- QC_Plots_Mito(seurat_object = gbm, high_cutoff = 12, pt.size = 0)
wrap_plots(p1, p2, p3, ncol = 3)
dev.off()

# remove low quality cells
gbm = subset(gbm, subset = nCount_RNA > 1000 & nCount_RNA < 100000 &
               nFeature_RNA > 500 & nFeature_RNA < 10000 &
               percent_mito < 12)

# Pass the features of interest as the `vars` parameter to FetchData
scatter_data = FetchData(gbm, vars = c('percent_GRCm39', 'percent_GRCh38',
                                       'GRCh38_exp', 'GRCm39_exp', 'Luciferase'))

p = ggplot(data = scatter_data) +
  geom_point(mapping = aes(x = GRCm39_exp, y = GRCh38_exp, 
                           color = Luciferase), 
             size = 0.1) +
  scale_color_gradientn(colors = c('orange', 'blue')) +
  cowplot::theme_cowplot()

ggsave(filename = 'scatter_GRCm39_exp_GRCh38_exp_color_Luciferase1.pdf', 
       plot = p, device = 'pdf', width = 4.5, height = 4, units = 'in')

pdf(file = 'vlnplot_human_vs_mouse_cells2.pdf', width = 5, height = 6)
VlnPlot(gbm, features = c('percent_GRCm39', 'percent_GRCh38', 'Luciferase'), 
        raster = F, ncol = 3, pt.size = 0) 

VlnPlot(gbm, features = c('GRCm39_exp', 'GRCh38_exp', 'Luciferase'),
        raster = F, ncol = 3, pt.size = 0)
dev.off()

pdf(file = 'vlnplot_percent_GRCh38_1.pdf', width = 4, height = 4)
VlnPlot(gbm, features = c('percent_GRCh38'), raster = F, pt.size = 0) +
  geom_hline(yintercept = c(10, 90), linetype = "dashed", color = "red")
dev.off()

# assign human and mouse cells ----
human_cells = subset(gbm, subset = Luciferase > 0 | percent_GRCh38 > 90)
mouse_cells = subset(gbm, subset = Luciferase == 0 & percent_GRCh38 < 10)

pdf(file = 'vlnplot_luciferase_expressed_human_aligned_percent1.pdf', 
    width = 6, height = 6)
pdf(file = 'vlnplot_qc_mouse_cells1.pdf', 
    width = 6, height = 6)
VlnPlot(human_cells, features = c('percent_GRCm39', 'percent_GRCh38', 'Luciferase'))
VlnPlot(mouse_cells, features = c('percent_GRCm39', 'percent_GRCh38', 'Luciferase'))
dev.off()

human_cells@active.ident = human_cells$orig.ident
mouse_cells@active.ident = mouse_cells$orig.ident

pdf(file = 'qc_vlnPlot_human_cells_1_1.pdf', width = 6, height = 4)
VlnPlot(human_cells, 
        features = c('nCount_RNA', 'nFeature_RNA', 'percent_mito'), 
        raster = F, ncol = 3, pt.size = 0)
dev.off()

saveRDS(human_cells, 'human_cells.rds')
saveRDS(mouse_cells, 'mouse_cells.rds')
saveRDS(gbm, 'gbm_filtered_qc.rds')

gbm = readRDS('gbm_filtered_qc.rds')
gbm[['Luc_pos']] = ifelse(GetAssayData(gbm, slot = 'counts')['Luciferase', ] > 0, 
                          'Luc+', 'Luc-')
Luc_pos = subset(gbm, subset = Luc_pos == 'Luc+')
Luc_pos[['human_mouse_refs']] = ifelse(Luc_pos$percent_GRCh38 < 10, 'mouse', 
                                       ifelse(Luc_pos$percent_GRCh38 > 90, 'human',
                                              NA))
pdf(file = 'pie_luciferase_pos_human_or_mouse.pdf', width = 5, height = 5)
pie(table(Luc_pos$human_mouse_refs))
dev.off()

## merge human and mouse cells passed qc ----
combined = merge(human_cells, y = mouse_cells,
                 project = 'GBM12_GBM38_HF2303_human_mouse_cells')
## Batch Correction ------------------------------------------------------------
# Following harmony seurat vignettes
# https://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/doc/Seurat.html
combined <- combined %>%
  NormalizeData(verbose = FALSE)
# Calculate variable features for each sample separately
VariableFeatures(combined) <- split(row.names(combined@meta.data), 
                                    combined@meta.data$orig.ident) %>% 
  lapply(function(cells_use) {
    combined[, cells_use] %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
      VariableFeatures()
  }) %>% unlist %>% unique

# performs scaling on the previously identified variable features (defaults to 2000)
# scaled data is used as input to PCA
combined <- combined %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = VariableFeatures(combined), npcs = 35, verbose = FALSE)

ElbowPlot(combined, ndims = 35, reduction = 'pca')

# batch correction on samples
combined = RunHarmony(combined, c("orig.ident"), max.iter.harmony = 50, 
                      plot_convergence = TRUE)
# harmony converged after 9 iterations

# Inspection of the modalities 
p1 = DimPlot_scCustom(combined, reduction = "harmony", pt.size = 0.1, group.by = "orig.ident")
p2 = DimPlot_scCustom(combined, reduction = "pca", pt.size = 0.1, group.by = "orig.ident")
pdf(file = 'dimplot_qc_harmony.pdf', width = 14, height = 4)
cowplot::plot_grid(p2, p1)
dev.off()

p2 = VlnPlot_scCustom(combined, features = "harmony_1", group.by = "orig.ident",  pt.size = 0)
p1 = VlnPlot_scCustom(combined, features = "PC_1", group.by = "orig.ident",  pt.size = 0)
pdf(file = 'vlnplot_qc_harmony.pdf', width = 12, height = 4)
cowplot::plot_grid(p1,p2)
dev.off()

p2 = VlnPlot_scCustom(combined, features = "harmony_2", group.by = "orig.ident",  pt.size = 0)
p1 = VlnPlot_scCustom(combined, features = "PC_2", group.by = "orig.ident",  pt.size = 0)
pdf(file = 'vlnplot_qc_harmony2.pdf', width = 12, height = 4)
cowplot::plot_grid(p1,p2)
dev.off()

# UMAP
# The default method for RunUMAP has changed from calling Python UMAP 
# via reticulate to the R-native UWOT using the cosine metric
combined <- combined %>%
  FindNeighbors(reduction = "harmony", dims = 1:35) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(reduction = "harmony", dims = 1:35)

gc()

save_figure <- function(reduction, group, object = combined, 
                        object_name = 'combined'){
  # save a 2D visualization of gbm object for a column of metadata
  pdf(file = paste(object_name, reduction, group, '2.pdf', sep = '_'), 
      width = 8, height = 6)
  p <- DimPlot_scCustom(object, reduction = reduction, 
                        group.by = group, label = F, 
                        ggplot_default_colors = F,
                        raster = F, pt.size = 0.1)
  print(p)
  dev.off()
}

save_figure('umap', 'seurat_clusters')
save_figure('umap', 'orig.ident')

markers <- FindAllMarkers(combined, only.pos = TRUE)
saveRDS(combined, 'combined_res01.rds')
saveRDS(markers, file = 'markers_res_01.rds')

markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  print(n = 200) -> top20

pdf(file = 'heatmap_top20.pdf', width = 9, height = 20)
DoHeatmap(combined, features = top20$gene) + 
  scale_fill_gradientn(colours = rev(brewer.pal(11, name = 'PiYG'))) 
dev.off()

pdf(file = 'featureplt_qc_combined_mito.pdf', width = 15, height = 5)
FeaturePlot(combined, features = c('nCount_RNA', 'nFeature_RNA', 'percent_mito'), 
            ncol = 3)
dev.off()

# a more refined resolution
combined = FindClusters(combined, resolution = 0.5)

pdf(file = 'dimplot_res05.pdf', width = 8, height = 6)
DimPlot_scCustom(combined, reduction = 'umap')
dev.off()

human_cell_ids = human_cells@assays$RNA@data@Dimnames[[2]]
mouse_cell_ids = mouse_cells@assays$RNA@data@Dimnames[[2]]

combined[['human_vs_mouse']] = NA
combined@meta.data[human_cell_ids, 'human_vs_mouse'] = 'human'
combined@meta.data[mouse_cell_ids, 'human_vs_mouse'] = 'mouse'
save_figure('umap', 'human_vs_mouse')

combined[['brain_vs_tumor']] = NA
combined@meta.data[combined$orig.ident == 'GBM12_brain_12000_SP2' |
                     combined$orig.ident == 'GBM38_brain_11913_SP2' |
                     combined$orig.ident == 'HF2303_brain_12002_SP2', 
                   'brain_vs_tumor'] = 'brain'

combined@meta.data[combined$orig.ident == 'GBM12_tumor_12000_SP1' |
                     combined$orig.ident == 'GBM38_tumor_11913_SP1' |
                     combined$orig.ident == 'HF2303_tumor_12002_SP1', 
                   'brain_vs_tumor'] = 'tumor'
save_figure('umap', 'brain_vs_tumor')

combined[['sample_origin_h_m']] = paste(combined$human_vs_mouse, 
                                        combined$brain_vs_tumor, sep = '_')
  
save_figure('umap', 'sample_origin_h_m')


combined[['res.0.1_refined_cl_5']] = as.character(combined$RNA_snn_res.0.1)
combined@meta.data[combined$RNA_snn_res.0.1 == 5 & 
                     (combined$sample_origin_h_m == 'mouse_brain' |
                        combined$sample_origin_h_m == 'mouse_tumor'), 
                   'res.0.1_refined_cl_5'] = "10"
combined$res.0.1_refined_cl_5 = factor(combined$res.0.1_refined_cl_5, 
                                       levels = c('0', '1', '2', '3', '4', '5', 
                                                  '6', '7', '8', '9', '10'))
save_figure('umap', 'res.0.1_refined_cl_5')

combined@active.ident = combined$res.0.1_refined_cl_5
levels(combined)
markers_res01_cl5 = FindAllMarkers(combined, only.pos = T)

saveRDS(markers_res01_cl5, file = 'markers_res01_cl5.rds')

markers_res01_cl5 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  print(n = 220) -> top20_cl5

pdf(file = 'heatmap_top20_cl5.pdf', width = 9, height = 21)
DoHeatmap(combined, features = top20_cl5$gene) + 
  scale_fill_gradientn(colours = rev(brewer.pal(11, name = 'PiYG'))) 
dev.off()


pdf(file = 'combined_dotplot_markers_cell_anno.pdf', width = 6, height = 11)
DotPlot_scCustom(combined, features = c(
  'Luciferase',
  'GRCh38-EGFR', 'GRCh38-GFAP', 'GRCh38-AQP4', 'GRCh38-PDGFRA', 'GRCh38-BCAN',
  'GRCh38-VCAN', 'GRCh38-GAD1', 'GRCh38-CDK4', # Neoplastic
  'GRCh38-SOX2', 'GRCh38-SOX4', 'GRCh38-TOP2A', 'GRCh38-MKI67', # proliferative
  'GRCm39-Slc1a2', 'GRCm39-Clu', 'GRCm39-Gfap', 'GRCm39-Cpe', 'GRCm39-Slc1a3', 
  'GRCm39-Aldoc', 'GRCm39-Aldh1l1', 'GRCm39-Ntsr2',# Astrocyte
  'GRCm39-Reln', 'GRCm39-Gad2', 'GRCm39-Gad1', 'GRCm39-Snap25', 'GRCm39-Bcl11b', 
  'GRCm39-Grin1', # Neuron
  'GRCm39-Cx3cr1', 'GRCm39-Tyrobp', 'GRCm39-Tgfbi', 'GRCm39-Ctss', 'GRCm39-Ptprc',
  'GRCm39-P2ry12', # Myeloid
  'GRCm39-Pecam1', 'GRCm39-Apold1', 'GRCm39-Tek', # Vascular
  'GRCm39-Mog', 'GRCm39-Cldn11', 'GRCm39-Plp1', 'GRCm39-Opalin', # Oligodendrocyte
  'GRCm39-Vcan', 'GRCm39-Pdgfra', 'GRCm39-Gpr17', 'GRCm39-Cntn1', # OPC
  'GRCm39-Cd3d', 'GRCm39-Mzb1', 'GRCm39-Ms4a1', # Lymphoid
  'GRCm39-Nkg7', 'GRCm39-Ncr1', 'GRCm39-Klrd1' # NK
), flip_axes = T, x_lab_rotate = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))
dev.off()

## Refine clusters 2 and 10 ----------------------------------------------------
cl210 = subset(combined, idents = c('2', '10'))

# batch correction on samples
cl210 = RunHarmony(cl210, c("orig.ident"), max.iter.harmony = 50, 
                   plot_convergence = TRUE)
cl210 <- cl210 %>%
  FindNeighbors(reduction = "harmony", dims = 1:35) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(reduction = "harmony", dims = 1:35)

save_figure('umap', 'seurat_clusters', object = cl210, 
            object_name = 'cl210')
save_figure('umap', 'orig.ident', object = cl210, object_name = 'cl210')


pdf(file = 'cl210_dotplot_markers_cell_anno_1.pdf', width = 5, height = 7)
DotPlot_scCustom(cl210, features = c(
  'GRCm39-Slc1a2', 'GRCm39-Gfap', 'GRCm39-Aqp4', 'GRCm39-Clu', 'GRCm39-Cpe', 
  'GRCm39-Slc1a3', 
  'GRCm39-Aldoc', 'GRCm39-Aldh1l1', 'GRCm39-Ntsr2',# Astrocyte
  'GRCm39-Reln', 'GRCm39-Gad2', 'GRCm39-Gad1', 'GRCm39-Bcl11b', 
  'GRCm39-Grin1', 'GRCm39-Fstl5', 'GRCm39-Synpr', 'GRCm39-Gabrg2', 
  'GRCm39-Rbfox3', 'GRCm39-Slc6a17', 
  'GRCm39-Cdh9', 'GRCm39-Kcnc2', 
  'GRCm39-Frmpd4', # Neuron
  'GRCm39-Mog', 'GRCm39-Cldn11', 'GRCm39-Plp1', 'GRCm39-Opalin', # Oligodendrocyte
  'GRCm39-Vcan', 'GRCm39-Pdgfra', 'GRCm39-Gpr17' # OPC
), flip_axes = T, x_lab_rotate = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')), 
                         limits = c(-2.5, 2.5))
dev.off()


cl210 = FindClusters(cl210, resolution = 1.2)
save_figure('umap', 'seurat_clusters', object = cl210, object_name = 'cl210')

cl210 = SetIdent(cl210, value = 'RNA_snn_res.0.8')
levels(cl210)

# expression of marker genes in combined clusters2&10 (now 13 clusters)
pdf(file = 'cl210_dotplot_markers_cell_anno.pdf', width = 7, height = 8)
DotPlot_scCustom(cl210, features = paste('GRCm39-', c(
  'Slc1a2', 'Aqp4', 'Gja1', 'Slco1c1', 'Aldh1l1', 'Slc1a3', 'Clu', # Astrocyte
  'Slc17a7', 'Enc1', 'Chn1', 'Pcp4', 'Gad2', 'Gad1',  'Reln', 'Tafa2', 
  'Ntng1', 'Fstl4', 'Grm1', 'Eml5', 'Camk4', # Neuron
  'Mog', 'Mag', 'Cldn11', 'Plp1', # Oligodendrocyte
  'Olig1', 'Vcan', 'Pdgfra', 'Gpr17', #OPC
  'Cfap54', 'Dnah9', 'Cfap299', 'Ttr', 'Prlr', # Ependymal
  'Cx3cr1', 'Tyrobp', 'Tgfbi', 'Ctss', 'Ptprc', 'P2ry12' # Myeloid
  ), sep = ''), 
  flip_axes = T, x_lab_rotate = T, group.by = 'RNA_snn_res.0.8') +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')))
dev.off()

markers_cl210 = FindAllMarkers(cl210, only.pos = T)
saveRDS(markers_cl210, 'markers_cl210_res08.rds')
markers_cl210 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  print(n = 280) -> top20cl210

pdf(file = 'heatmap_top20_cl210.pdf', width = 9, height = 21)
DoHeatmap(cl210, features = top20cl210$gene) + 
  scale_fill_gradientn(colours = rev(brewer.pal(11, name = 'PiYG'))) 
dev.off()

# cell annotations of clusters2&10 
cl210 = RenameIdents(cl210, 
                     '3' = 'Astrocyte', '7' = 'Astrocyte', 
                     '11' = 'Astrocyte', 
                     '0' = 'Neuron_1', '1' = 'Neuron_1', 
                     '2' = 'Neuron_1', '4' = 'Neuron_1', 
                     '5' = 'Neuron_1', '6' = 'Neuron_1', 
                     '8' = 'Neuron_2', '9' = 'Neuron_2', 
                     '10' = 'Oligodendrocyte', '12' = 'OPC', 
                     '13' = 'Ependymal/Myeloid')
cl210[['cell_anno']] = cl210@active.ident

# cluster 13 has mixture of ependymal and myeloid signature
cl13 = subset(cl210, idents = '13')
cl13 <- cl13 %>%
  FindNeighbors(reduction = "harmony", dims = 1:35) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(reduction = "harmony", dims = 1:35)

pdf(file = 'cl13_dotplot_markers_cell_anno.pdf', width = 5, height = 4)
DotPlot_scCustom(cl13, features = paste('GRCm39-', c('Cx3cr1', 'Tyrobp', 'Tgfbi', 
                                       'Ctss', 'Ptprc', 'P2ry12', # Myeloid
                                       'Cfap54', 'Dnah9', 'Cfap299', 
                                       'Ttr', 'Prlr'), sep = ''), # Ependymal
                 flip_axes = T, x_lab_rotate = T, group.by = 'RNA_snn_res.0.8') +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')))
dev.off()
# cluster 1 --> Myeloid, cluster 2 --> Ependymal, 0: Not Determined
cl13 = RenameIdents(cl13, '0' = 'Not Determined', 
                    '1' = 'Myeloid', '2' = 'Ependymal')
cl13[['cell_anno']] = cl13@active.ident
saveRDS(cl13, 'cl13.rds')

# add cluster 13 annotation to clusters2&10
cellIds_cl13 = colnames(cl13)

cl210[['cell_anno']] = as.character(cl210$cell_anno)
cl210@meta.data[cellIds_cl13, "cell_anno"] = as.character(cl13$cell_anno)
cl210$cell_anno = factor(cl210$cell_anno,
                         levels = c('Neuron_1', 'Neuron_2', 'Astrocyte', 
                                    'Oligodendrocyte', 'OPC', 
                                    'Myeloid', 'Ependymal', 'Not Determined'))

save_figure('umap', 'cell_anno', object = cl210, object_name = 'cl210')

pdf(file = 'cl210_dotplot_serine_genes.pdf', width = 7, height = 8)
DotPlot_scCustom(cl210, features = serine_genes_mouse, 
                 flip_axes = T, x_lab_rotate = T, group.by = 'cell_anno') +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')))

dev.off()

saveRDS(cl210, 'cl210.rds')

## add cell annotations of clusters 2 and 10 to the combined object ------------
combined = readRDS('../combined_res01.rds')

combined = SetIdent(combined, value = 'res.0.1_refined_cl_5')
combined = RenameIdents(combined, '0' = 'Neoplastic', '1' = 'Neoplastic', 
                        '3' = 'Neoplastic', '5' = 'Neoplastic', '8' = 'Neoplastic',
                        '4' = 'Myeloid', '7' = 'Myeloid', 
                        '6' = 'Vascular', '9' = 'Oligodendrocyte',
                        '2' = 'Glial/Neuronal', '10' = 'Glial/Neuronal')

combined[['coarse_cell_anno']] = as.character(combined@active.ident)
combined$coarse_cell_anno = factor(combined$coarse_cell_anno,
                                   levels = c('Neoplastic', 'Glial/Neuronal',
                                              'Oligodendrocyte', 'Myeloid',
                                              'Vascular'))
save_figure('umap', 'coarse_cell_anno', object = combined, 
            object_name = 'combined')

cellIds_cl210 = colnames(cl210)

combined[['coarse_cell_anno_2']] = as.character(combined$coarse_cell_anno)
combined@meta.data[cellIds_cl210, 
                   "coarse_cell_anno_2"] = as.character(cl210$cell_anno)
levels(factor(combined$coarse_cell_anno_2))
combined$coarse_cell_anno_2 = factor(combined$coarse_cell_anno_2, 
                                     levels = c('Neoplastic', 
                                                'Neuron_1', 'Neuron_2',
                                                'Astrocyte', 
                                                'Oligodendrocyte', 
                                                'OPC', 
                                                'Myeloid', 
                                                'Ependymal',
                                                'Vascular',
                                                'Not Determined'))

save_figure('umap', 'coarse_cell_anno_2', object = combined, 
            object_name = 'combined')
combined[['cell_anno']] = as.character(combined$coarse_cell_anno_2)
combined@meta.data[(combined$cell_anno == 'Neuron_1') | 
                     (combined$cell_anno == 'Neuron_2'),
                   'cell_anno'] = 'Neuron'
combined$cell_anno = factor(combined$cell_anno,
                            levels = c('Neoplastic', 'Neuron', 'Astrocyte', 
                                       'Oligodendrocyte', 'OPC', 'Myeloid', 
                                       'Ependymal', 'Vascular',
                                       'Not Determined'))
save_figure('umap', 'cell_anno', object = combined, object_name = 'combined')

saveRDS(combined, 'combined_res01.rds')

data("palette36")
# Proportion plot
plot_integrated_clusters = function (srat, cluster_x, cluster_y, is_cluster_x_int, 
                                     filename) { 
  # take an integrated Seurat object, plot distributions over orig.ident
  
  count_table <- table(srat@meta.data[[cluster_x]], srat@meta.data[[cluster_y]])
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)
  
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
    scale_fill_manual(values = unname(palette36)) +
    ylab("Fraction of cells in each dataset") + xlab("Cluster number") + 
    theme(legend.position="top")
  
  p = p2 + p1 + plot_layout(widths = c(3, 1))
  ggsave(filename = paste(filename, '.pdf', sep = ''),
         plot = p, device = 'pdf', width = 12, height = 5.5, units = 'in')
  
}

plot_integrated_clusters(combined, 
                         cluster_x = 'human_vs_mouse', 
                         cluster_y = 'coarse_cell_anno_2', 
                         is_cluster_x_int = F, 
                         filename = 'human_vs_mouse_coarse_cell_anno_2')
combined$human_vs_mouse

combined[['level_1_anno']] = as.character(combined$coarse_cell_anno_2)

levels(combined$coarse_cell_anno_2)
levels(factor(combined$level_1_anno))
combined@meta.data[combined$coarse_cell_anno_2 != 'Neoplastic', 
                   'level_1_anno'] = 'Non-Neoplastic'
combined[['human_vs_mouse_l1']] = paste(combined$human_vs_mouse, 
                                        combined$level_1_anno, sep = '_')
save_figure('umap', 'human_vs_mouse_l1', object = combined, 
            object_name = 'combined')

table(combined$human_vs_mouse_l1)
# human_Neoplastic human_Non-Neoplastic mouse_Non-Neoplastic 
# 12829                  158                 7070
saveRDS(combined, 'combined_res01_updatedcl210.rds')

# human cells are all neoplastic ----
setwd('./human_neoplastic_mouse_non_neoplastic')
combined@active.ident = factor(combined$human_vs_mouse_l1)
combined = subset(combined, idents = 'human_Non-Neoplastic', invert = T)


combined[['mouse_model']] = as.character(combined$orig.ident)
combined@meta.data[(combined$orig.ident == 'HF2303_tumor_12002_SP1') |
                     (combined$orig.ident == 'HF2303_brain_12002_SP2'),
                   'mouse_model'] = 'HF2303'
combined@meta.data[(combined$orig.ident == 'GBM38_tumor_11913_SP1') |
                     (combined$orig.ident == 'GBM38_brain_11913_SP2'),
                   'mouse_model'] = 'GBM38'
combined@meta.data[(combined$orig.ident == 'GBM12_tumor_12000_SP1') |
                     (combined$orig.ident == 'GBM12_brain_12000_SP2'),
                   'mouse_model'] = 'GBM12'

# update plots for filtered data

save_figure('umap', 'res.0.1_refined_cl_5', object = combined, 
            object_name = 'combined_filtered_')

plot_integrated_clusters(combined, 
                         cluster_x = 'coarse_cell_anno_2',
                         cluster_y = 'mouse_model', 
                         is_cluster_x_int = F, 
                         filename = 'coarse_cell_anno_2_mouse_model')

saveRDS(combined, 'combined.rds')
with(combined@meta.data, table(orig.ident, cell_anno))

# multiple line comment: ctrl + shift + C

# orig.ident               Neoplastic Neuron_1 Neuron_2 Astrocyte Oligodendrocyte  OPC
# GBM12_brain_12000_SP2         105     1961      138       406              99   19
# GBM12_tumor_12000_SP1        2992        5        3        24              19    9
# GBM38_brain_11913_SP2          12      639       73       253              90   33
# GBM38_tumor_11913_SP1        5950        0        1         1               3   12
# HF2303_brain_12002_SP2         64      386       61       180              27    2
# HF2303_tumor_12002_SP1       3706        0        3         5              16    3

# orig.ident               Myeloid Ependymal Vascular
# GBM12_brain_12000_SP2       94        12      101
# GBM12_tumor_12000_SP1      268         7       72
# GBM38_brain_11913_SP2      317        14       93
# GBM38_tumor_11913_SP1      817         9       32
# HF2303_brain_12002_SP2     133        13      108
# HF2303_tumor_12002_SP1     414         7       88

# coarse_cell_anno_2
# orig.ident               Neoplastic Neuron_1 Neuron_2 Astrocyte Oligodendrocyte  OPC
# GBM12_brain_12000_SP2         105     1961      138       406              99   19
# GBM12_tumor_12000_SP1        2992        5        3        24              19    9
# GBM38_brain_11913_SP2          12      639       73       253              90   33
# GBM38_tumor_11913_SP1        5950        0        1         1               3   12
# HF2303_brain_12002_SP2         64      386       61       180              27    2
# HF2303_tumor_12002_SP1       3706        0        3         5              16    3
# coarse_cell_anno_2
# orig.ident               Myeloid Ependymal Vascular Not Determined
# GBM12_brain_12000_SP2       94         4      101              8
# GBM12_tumor_12000_SP1      271         2       72              2
# GBM38_brain_11913_SP2      317         3       93             11
# GBM38_tumor_11913_SP1      826         0       32              0
# HF2303_brain_12002_SP2     133         2      108             11
# HF2303_tumor_12002_SP1     414         3       88              4

# cell_anno
# orig.ident               Neoplastic Neuron Astrocyte Oligodendrocyte  OPC Myeloid
# GBM12_brain_12000_SP2         105   2099       406              99   19      94
# GBM12_tumor_12000_SP1        2992      8        24              19    9     271
# GBM38_brain_11913_SP2          12    712       253              90   33     317
# GBM38_tumor_11913_SP1        5950      1         1               3   12     826
# HF2303_brain_12002_SP2         64    447       180              27    2     133
# HF2303_tumor_12002_SP1       3706      3         5              16    3     414
# cell_anno
# orig.ident               Ependymal Vascular Not Determined
# GBM12_brain_12000_SP2          4      101              8
# GBM12_tumor_12000_SP1          2       72              2
# GBM38_brain_11913_SP2          3       93             11
# GBM38_tumor_11913_SP1          0       32              0
# HF2303_brain_12002_SP2         2      108             11
# HF2303_tumor_12002_SP1         3       88              4
