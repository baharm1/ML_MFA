# Preprocess our scRNA-seq data

# Import Libraries ------------------------------------------------------------
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scCustomize)
library(stringr)

# Directories -----------------------------------------------------------------
directory = './Baharan Meghdadi/scRNAseq_data_TRP_GBM38/11913-SP/cellranger_BM'
fig_directory = './Baharan Meghdadi/Baharan-Deepak/ML_manuscript/codes/single_cell_analysis/TRP_GBM38_files'

# Preprocessing Steps and Quality Control -------------------------------------

# Sample 3 (TRP-tumor) includes mouse cells -----------------------------------
setwd(directory)
sample_dir = list.files(path = '.', 
                        pattern = 'SP', 
                        all.files = T, full.names = T)
sample_names = unlist(strsplit(sample_dir, '/'))[c(FALSE, TRUE)]

sample_no = 3

barcodes_dir = list.files(path = paste(sample_dir[sample_no], 
                                       '/filtered_feature_bc_matrix/', 
                                       sep = ''),
                          pattern = 'barcodes.tsv.gz', 
                          all.files = T,
                          full.names = T)
barcodes = read.table(gzfile(barcodes_dir), stringsAsFactors = F)

matrix_dir = list.files(path = paste(sample_dir[sample_no], 
                                     '/filtered_feature_bc_matrix/', 
                                     sep = ''),
                        pattern = 'matrix.mtx.gz', 
                        all.files = T,
                        full.names = T)
mtrx = readMM(gzfile(matrix_dir))

features_dir = list.files(path = paste(sample_dir[sample_no], 
                                       '/filtered_feature_bc_matrix/', 
                                       sep = ''),
                          pattern = 'features.tsv.gz', 
                          all.files = T,
                          full.names = T)
features = read.table(gzfile(features_dir), stringsAsFactors = F)

dimnames(mtrx) = list(features$V2, 
                      paste(sample_names[sample_no], barcodes$V1, sep = '-'))
mtrx = CreateSeuratObject(counts = mtrx, project = sample_names[sample_no])

mtrx[['percent.mito']] = PercentageFeatureSet(mtrx, pattern = "^mt-")


pdf(file = 'qc_vlnPlot_before_SP3.pdf', width = 6, height = 4)
VlnPlot(mtrx, 
        features = c('nCount_RNA', 'nFeature_RNA', 'percent.mito'), 
        raster = F, ncol = 3, pt.size = 0)
dev.off()

setwd(fig_directory)

# remove low quality cells
mtrx = subset(mtrx, subset = nCount_RNA > 1000 & nCount_RNA < 100000 &
               nFeature_RNA > 500 & nFeature_RNA < 10000 &
               percent.mito < 12)

mtrx <- mtrx %>%
  NormalizeData(verbose = FALSE)
# Calculate variable features for each sample separately
VariableFeatures(mtrx) <- split(row.names(mtrx@meta.data), 
                                mtrx@meta.data$orig.ident) %>% 
  lapply(function(cells_use) {
    mtrx[, cells_use] %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
      VariableFeatures()
  }) %>% unlist %>% unique

# performs scaling on the previously identified variable features (defaults to 2000)
# scaled data is used as input to PCA
mtrx <- mtrx %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = VariableFeatures(mtrx), npcs = 35, verbose = FALSE)

ElbowPlot(mtrx, ndims = 35, reduction = 'pca')

# UMAP
# The default method for RunUMAP has changed from calling Python UMAP 
# via reticulate to the R-native UWOT using the cosine metric
mtrx <- mtrx %>%
  FindNeighbors(reduction = "pca", dims = 1:35) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(reduction = "pca", dims = 1:35)

gc()

save_figure <- function(reduction, group, object = mtrx, 
                        object_name = 'TRP'){
  # save a 2D visualization of gbm object for a column of metadata
  pdf(file = paste(object_name, reduction, group, '2507.pdf', sep = '_'), 
      width = 8, height = 6)
  p <- DimPlot_scCustom(object, reduction = reduction, 
                        group.by = group, label = F, 
                        ggplot_default_colors = F,
                        raster = F, pt.size = 0.1)
  print(p)
  dev.off()
}

save_figure('umap', 'seurat_clusters')
FeaturePlot(mtrx, features = c("percent.mito"), reduction = 'umap')
FeaturePlot(mtrx, features = c("Luciferase"), reduction = 'umap')

markers <- FindAllMarkers(mtrx, only.pos = TRUE)

saveRDS(mtrx, 'SP3_res01_2507.rds')
saveRDS(markers, file = 'markers_SP3_res_01_2507.rds')

markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  print(n = 260) -> top20

mtrx = FindClusters(mtrx, resolution = 0.5)
save_figure('umap', 'seurat_clusters')

pdf(file = 'SP3_res05_Luciferase_2507.pdf', height = 10, width = 4)
DotPlot(mtrx, features = 'Luciferase')
dev.off()

pdf(file = 'SP3_dotplot_markers_cell_anno_250729.pdf', width = 6, height = 6)
DotPlot_scCustom(mtrx, features = c(
  'Luciferase',
  'Sox2', 'Sox4', 'Top2a', 'Mki67', # proliferative
  'Slc1a3', 'Aldoc', 
  'Cx3cr1', 'Tyrobp', 'Tgfbi', 'Ctss', 'P2ry12', # Myeloid
  'Pecam1', 'Apold1', 'Tek', # Vascular
  'Mog', 'Cldn11', 'Plp1', 'Opalin', # Oligodendrocyte
  'Vcan', 'Pdgfra', # OPC
  'Cd3d', 'Nkg7', 'Ncr1', 'Klrd1', 'Mzb1', 'Ms4a1' # Lymphoid
), flip_axes = T, x_lab_rotate = T) + #group.by = 'RNA_snn_res.0.1', 
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')),
                         limits = c(-2.5, 2.5))
dev.off()

cl2 = subset(mtrx, idents = '2')
cl2 <- cl2 %>%
  FindNeighbors(reduction = "pca", dims = 1:35) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(reduction = "pca", dims = 1:35)

save_figure('umap', 'seurat_clusters', object = cl2, 
            object_name = 'SP3_cl2')

pdf(file = 'SP3_cl2_dotplot_markers_cell_anno_2507_2.pdf', width = 7, height = 8)
DotPlot_scCustom(cl2, features = c(
  'Slc1a2', 'Aqp4', 'Gja1', 'Slco1c1', 'Aldh1l1', 'Slc1a3', 'Clu', # Astrocyte
  'Slc17a7', 'Enc1', 'Chn1', 'Pcp4', 'Gad2', 'Gad1',  'Reln', 'Tafa2', 
  'Ntng1', 'Fstl4', 'Grm1', 'Eml5', 'Camk4', # Neuron
  'Mog', 'Mag', 'Cldn11', 'Plp1', # Oligodendrocyte
  'Olig1', 'Vcan', 'Pdgfra', 'Gpr17', #OPC
  'Cfap54', 'Dnah9', 'Cfap299', 'Ttr', 'Prlr', # Ependymal
  'Cx3cr1', 'Tyrobp', 'Tgfbi', 'Ctss', 'Ptprc', 'P2ry12' # Myeloid
), 
flip_axes = T, x_lab_rotate = T, group.by = 'RNA_snn_res.0.8') +
  scale_colour_gradientn(colours = rev(brewer.pal(11, name = 'RdYlBu')))
dev.off()

cl2 = SetIdent(cl2, value = 'RNA_snn_res.0.8')
# cell annotations of cluster 2
cl2 = RenameIdents(cl2, '0' = 'Neoplastic', '1' = 'Neoplastic', 
                   '3' = 'Neoplastic', '4' = 'Neoplastic', 
                   '5' = 'Neoplastic', '7' = 'Neoplastic', 
                   '2' = 'Neoplastic', '6' = 'Neoplastic', '9' = 'Neoplastic',
                   '8' = 'Neoplastic',  
                   '10' = 'Oligodendrocyte')
cl2[['cell_anno']] = cl2@active.ident

save_figure('umap', 'cell_anno', object = cl2, 
            object_name = 'SP3_cl2_2')

# cell annotations of SP3
mtrx = RenameIdents(mtrx, 
                    '0' = 'Myeloid', '1' = 'Lymphoid', '2' = 'Mixed',
                    '3' = 'Neoplastic', '4' = 'Myeloid',
                    '5' = 'Lymphoid', '6' = 'Lymphoid', 
                    '7' = 'Myeloid', '8' = 'Myeloid',
                    '9' = 'Myeloid', '10' = 'Myeloid',
                    '11' = 'Vascular', '12' = 'OPC')

mtrx[['cell_anno']] = as.character(mtrx@active.ident)

cellIds_cl2 = colnames(cl2)

mtrx[['cell_anno_2']] = mtrx$cell_anno

mtrx@meta.data[cellIds_cl2, 
                   "cell_anno_2"] = as.character(cl2$cell_anno)
levels(factor(mtrx$cell_anno_2))

mtrx$cell_anno_2 = factor(mtrx$cell_anno_2, 
                                     levels = c('Neoplastic', 
                                                'Oligodendrocyte', 
                                                'OPC', 
                                                'Myeloid', 
                                                'Lymphoid',
                                                'Vascular'))

save_figure('umap', 'cell_anno_2', object = mtrx, 
            object_name = 'SP3_2')

saveRDS(mtrx, 'SP3_res01_2507.rds')
