# Preprocess our scRNA-seq data

# Import Libraries ------------------------------------------------------------
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
library(Polychrome)
library(scCustomize)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(harmony)
library(GSEABase)

# Directories -----------------------------------------------------------------
directory = './baharan/reanalysis/'
fig_directory = './baharan/reanalysis/figs'
seurat_directory = './baharan/reanalysis/seu_obj'

# Plot Functions --------------------------------------------------------------
qc_visualization <- function(data, thr_umi, thr_gene, thr_mito, thr_ribo){
  # Visualize the number UMIs/transcripts per cell
  p = ggplot(data@meta.data, aes(x = nCount_RNA)) + 
    geom_density() + 
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = thr_umi)
  ggsave(filename = 'umi.pdf', 
         plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')
  
  # Visualize the distribution of genes detected per cell via histogram
  p = ggplot(data@meta.data, aes(x = nFeature_RNA)) + 
    geom_density() + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = thr_gene)
  ggsave(filename = 'gene.pdf', 
         plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')
  
  # Visualize the distribution of mitochondrial gene expression detected per cell
  p = ggplot(data@meta.data, aes(x = percent.mito)) + 
    geom_density() + 
    theme_classic() +
    geom_vline(xintercept = thr_mito)
  ggsave(filename = 'mito.pdf', 
         plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')
  
  # Visualize the distribution of mitochondrial gene expression detected per cell
  p = ggplot(data@meta.data, aes(x = percent.ribo)) + 
    geom_density() + 
    theme_classic() +
    geom_vline(xintercept = thr_ribo)
  ggsave(filename = 'ribo.pdf', 
         plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')
  
  # Visualize the correlation between genes detected and number of UMIs and 
  # determine whether strong presence of cells with low numbers of genes/UMIs
  p = ggplot(data@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mito)) + 
    geom_point(size = 0.1) + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_hline(yintercept = thr_gene) +
    geom_vline(xintercept = thr_umi)
  ggsave(filename = 'umi_gene.pdf', 
         plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')
}

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

save_figure <- function(reduction, group){
  # save a 2D visualization of gbm object for a column of metadata
  pdf(file = paste('wna', reduction, group, '.pdf', sep = '_'), 
      width = 8, height = 6)
  p <- DimPlot_scCustom(gbm, reduction = reduction, group.by = group, label = F, 
                        raster = F, pt.size = 0.1)
  print(p)
  dev.off()
}

# function of scrabble package
hierarchy = function(m, quadrants = NULL, log.scale = T) {
  
  if (!is.null(quadrants)) {
    stopifnot(all(unlist(quadrants) %in% colnames(m)))
    dat = as.data.frame(sapply(quadrants, function(col) do.call(pmax, list(as.data.frame(m[, col])))))
  } else {
    stopifnot(ncol(m) == 4)
    dat = as.data.frame(m)
  }
  
  rows = rownames(m)
  colnames(dat) = c('bl', 'br', 'tl', 'tr')
  
  dat = dat %>%
    dplyr::mutate(bottom = pmax(bl, br),
                  top = pmax(tl, tr),
                  b.center = br - bl,
                  t.center = tr - tl,
                  x = ifelse(bottom > top, b.center, t.center), # dependent var
                  x.scaled = (sign(x) * log2(abs(x) + 1)),
                  y = top - bottom, # independent var
                  y.scaled = (sign(y) * log2(abs(y) + 1)))
  
  if (!log.scale) dat = dplyr::transmute(dat, X = x, Y = y)
  else dat = dplyr::transmute(dat, X = x.scaled, Y = y.scaled)
  rownames(dat) = rows
  class(dat) = append(class(dat), 'hierarchy')
  dat
}

# function of scrabble package
plot_hierarchy = function(X,
                          quadrant.names = c('bl', 'br', 'tl', 'tr'),
                          main = NULL,
                          xlab = 'Relative meta-module score [log2(|SC1-SC2|+1)]',
                          ylab = 'Relative meta-module score [log2(|SC1-SC2|+1)]',
                          groups = NULL,
                          group.cols = NULL, 
                          legend = T,
                          legend.pos = 'bottom',
                          legend.horiz = T) {
  
  if (is.null(groups)) col = 'darkred'
      else col = 'grey85'
          
      cex = 0.3
      plot(X[,1], X[,2], pch = 20, cex = cex, col = col, main = main, 
           xlab = xlab, ylab = ylab)
      
      if (is.null(groups)) legend = F
      else {
        stopifnot(!is.null(names(groups)))
        stopifnot(all(groups %in% rownames(X)))
        groups = split(groups, names(groups))
        Xgrp = sapply(groups, function(rows) X[rows,,drop = F], simplify = F)
        if (!is.null(group.cols)) colgrp = group.cols[names(groups)]
        else colgrp = rainbow(n = length(Xgrp))
        Map(points,
            x = sapply(Xgrp, `[[`, 1, simplify = F),
            y = sapply(Xgrp, `[[`, 2, simplify = F),
            col = colgrp,
            MoreArgs = list(pch = 20, cex = cex))
      }
      
      abline(v = 0, lty = 2)
      abline(h = 0, lty = 2)
      
      if (legend) {
        legend(legend.pos,
               fill = colgrp,
               legend = names(groups),
               horiz = legend.horiz,
               cex = 0.8,
               box.col = 'white',
               bg = 'white',
               box.lwd = 0)
      }
      
      Names = quadrant.names
      cex = 1.2
      mtext(side = 1, adj = 0, text = Names[1], cex = cex, line = cex - 1)
      mtext(side = 1, adj = 1, text = Names[2], cex = cex, line = cex - 1)
      mtext(side = 3, adj = 0, text = Names[3], cex = cex)
      mtext(side = 3, adj = 1, text = Names[4], cex = cex)
}

# Preprocessing Steps and Quality Control -------------------------------------

## read filtered count matrices ----
setwd(directory)
sample_dir = list.files(path = '.', 
                        pattern = 'Sample_', 
                        all.files = T, full.names = T)
sample_names = unlist(strsplit(sample_dir, '/'))[c(FALSE, TRUE)]

mtrx_list = c()

for (sample_no in 1:length(sample_dir)){
  barcodes_dir = list.files(path = paste(sample_dir[sample_no], 
                                         '/outs/filtered_feature_bc_matrix/', 
                                         sep = ''),
                       pattern = 'barcodes.tsv.gz', 
                       all.files = T,
                       full.names = T)
  barcodes = read.table(gzfile(barcodes_dir), stringsAsFactors = F)
  
  matrix_dir = list.files(path = paste(sample_dir[sample_no], 
                                       '/outs/filtered_feature_bc_matrix/', 
                                       sep = ''),
                          pattern = 'matrix.mtx.gz', 
                          all.files = T,
                          full.names = T)
  mtrx = readMM(gzfile(matrix_dir))
  
  features_dir = list.files(path = paste(sample_dir[sample_no], 
                                         '/outs/filtered_feature_bc_matrix/', 
                                         sep = ''),
                           pattern = 'features.tsv.gz', 
                           all.files = T,
                           full.names = T)
  features = read.table(gzfile(features_dir), stringsAsFactors = F)
  
  dimnames(mtrx) = list(features$V2, 
                        paste(sample_names[sample_no], barcodes$V1, sep = '_'))
  mtrx = CreateSeuratObject(counts = mtrx, project = sample_names[sample_no])
  mtrx$orig.ident = sample_names[sample_no]
  mtrx_list = append(mtrx, mtrx_list)
  rm(mtrx, barcodes, features)
}

## merge filtered counts into a seurat object ----
gbm = merge(mtrx_list[[1]], y = mtrx_list[2:length(sample_names)],
            project = 'sc_wna')

# patient samples
gbm$orig.ident = factor(gbm$orig.ident)
gbm@active.ident = gbm$orig.ident

## save merged filtered counts ----
getwd()
setwd(seurat_directory)
saveRDS(gbm, 'gbm_merged_seurat.rds')
rm(mtrx_list)

# add mitochondrial and ribosomal transcripts to metadata
gbm[['percent.mito']] = PercentageFeatureSet(gbm, pattern = "^MT-")
gbm[['percent.ribo']] = PercentageFeatureSet(gbm, 
                                             pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

gbm = SetIdent(gbm, value = gbm@project.name)

setwd(fig_directory)
pdf(file = 'qc_vlnPlot_before.pdf', width = 7, height = 4)
VlnPlot(gbm, 
        features = c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'percent.ribo'), 
        raster = F, ncol = 4, pt.size = 0)
dev.off()

qc_features = colnames(gbm@meta.data)
qc_features = qc_features[-c(1)] # remove orig.ident
qc_features
min_qc_features = sapply(qc_features, function(x) min(gbm@meta.data[, x]))
max_qc_features = sapply(qc_features, function(x) max(gbm@meta.data[, x]))

# define thresholds for number of UMI counts, number of detected genes, 
# percentages of mitochondrial and ribosomal transcripts
thr_umi = 1000
thr_gene = 500
thr_mito = 12
thr_ribo = 50

qc_visualization(gbm, thr_umi, thr_gene, thr_mito, thr_ribo) 

# remove low quality cells
gbm_filtered = subset(gbm, 
                      subset = nCount_RNA > thr_umi & nFeature_RNA > thr_gene & 
                        percent.mito < thr_mito & percent.ribo < thr_ribo)

pdf(file = 'qc_vlnPlot_after.pdf', width = 7, height = 4)
VlnPlot(gbm_filtered, 
        features = c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'percent.ribo'), 
        raster = F, ncol = 4, pt.size = 0)
dev.off()

## Doublet Finder --------------------------------------------------------------
table(gbm_filtered$orig.ident)
gbm_split = SplitObject(gbm_filtered, split.by = 'orig.ident')

# multiplet rate from
# https://assets.ctfassets.net/an68im79xiti/4tjk4KvXzTWgTs8f3tvUjq/2259891d68c53693e753e1b45e42de2d/CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf
df = data.frame('multiplet_rate' = c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6),
                'no_cells_recovered' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))

# calculate regression line coefficents (slope and intercept)
mltp_l = lm(formula = multiplet_rate ~ no_cells_recovered, data = df)
coefs = unname(coefficients(mltp_l))

# plot regression line 
pdf('multiplet_rate.pdf', width = 5, height = 5)
plot(df$no_cells_recovered, df$multiplet_rate)
abline(lm(multiplet_rate ~ no_cells_recovered, data = df))
dev.off()

# DoubletFinder workflow
for (i in 1:length(gbm_split)){
  print(names(gbm_split)[i])
  
  # Seurat workflow
  gbm_sample = NormalizeData(gbm_split[[i]])
  gbm_sample = FindVariableFeatures(gbm_sample)
  gbm_sample = ScaleData(gbm_sample)
  gbm_sample = RunPCA(gbm_sample)
  print(ElbowPlot(gbm_sample))
  
  # Find significant PCs from 
  # https://rpubs.com/olneykimberly/LPS_brain_neuron_recluster
  stdv = gbm_sample[["pca"]]@stdev
  # percent of variation of PCs
  percent_stdv = (stdv / sum(stdv)) * 100
  # cumulative percent of variation for each PC
  cumulative = cumsum(percent_stdv)
  # Determine which PC exhibits cumulative percent greater than 90% and
  # and % variation associated with the PC as less than 5
  co1 = which(cumulative > 90 & percent_stdv < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 = sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                      percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  # PCs covering the majority of the variation in the data
  min_pc = min(co1, co2)
  print(min_pc)
  
  # clustering and visualization
  gbm_sample = RunUMAP(gbm_sample, dims = 1:min_pc)
  gbm_sample = FindNeighbors(object = gbm_sample, dims = 1:min_pc)              
  gbm_sample = FindClusters(object = gbm_sample, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep_list = paramSweep_v3(gbm_sample, PCs = 1:min_pc, num.cores = 1, sct = F)
  sweep_stats = summarizeSweep(sweep_list, GT = F)
  bcmvn = find.pK(sweep_stats)
  
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  bcmvn_max = bcmvn[which.max(bcmvn$BCmetric),]
  optimal_pk = bcmvn_max$pK
  optimal_pk = as.numeric(levels(optimal_pk))[optimal_pk]
  print(optimal_pk)
  
  # Homotypic doublet proportion estimate
  homotypic_prop = modelHomotypic(gbm_sample@meta.data$seurat_clusters) 
  
  # expected number of doublets
  # regressed multiplet rate 
  mltp = (coefs[1] + coefs[2] * ncol(gbm_sample)) / 100 
  print(mltp)
  nExp_poi = round(mltp * ncol(gbm_sample)) 
  nExp_poi_adj = round(nExp_poi * (1 - homotypic_prop))
  
  # run DoubletFinder
  gbm_sample = doubletFinder_v3(seu = gbm_sample, 
                                PCs = 1:min_pc, 
                                pK = optimal_pk,
                                nExp = nExp_poi_adj)
  
  # subset singlets
  colnames(gbm_sample@meta.data)[ncol(gbm_sample@meta.data)] <- "doublet_finder"
  table(gbm_sample$doublet_finder)
  gbm_split[[i]] = subset(gbm_sample, doublet_finder == "Singlet")
  
}
# merge filtered singlets into a Seurat object
gbm_singlets <- merge(x = gbm_split[[1]], y = gbm_split[2:length(gbm_split)],
                      project = "wna_singlets")
# save filtered singlets
setwd(seurat_directory)
saveRDS(gbm_singlets, 'seu_obj_2310/singlets.rds')
# check DoubletFinder output
gbm_singlets = SetIdent(gbm_singlets, value = gbm_singlets@project.name)
pdf(file = 'qc_vlnPlot_after_doublet_filtered.pdf', width = 7, height = 4)
VlnPlot(gbm_singlets, 
        features = c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'percent.ribo'), 
        raster = F, ncol = 4, pt.size = 0)
dev.off()
# remove potential doublets manually
gbm_singlets = subset(gbm_singlets, 
                      subset = nCount_RNA < 100000 & nFeature_RNA < 10000)
# save final filtered singlets
saveRDS(gbm_singlets, 'seu_obj_2310/singlets_filtered.rds')

## Batch Correction ------------------------------------------------------------

setwd(fig_directory)
gbm = readRDS('seu_obj_2310/singlets_filtered.rds')

gbm = SetIdent(gbm, value = gbm@project.name)

setwd(fig_directory)
pdf(file = 'qc_vlnPlot_singlets.pdf', width = 7, height = 4)
VlnPlot(gbm, 
        features = c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'percent.ribo'), 
        raster = F, ncol = 4, pt.size = 0)
dev.off()

qc_visualization(gbm, 0, 0, 0, 0)

gc()
# Following harmony seurat vignettes
# https://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/doc/Seurat.html
gbm <- gbm %>%
  NormalizeData(verbose = FALSE)
# Calculate variable features for each sample separately
VariableFeatures(gbm) <- split(row.names(gbm@meta.data), 
                               gbm@meta.data$orig.ident) %>% 
  lapply(function(cells_use) {
    gbm[, cells_use] %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    VariableFeatures()
    }) %>% unlist %>% unique

# performs scaling on the previously identified variable features (defaults to 2000)
# scaled data is used as input to PCA
gbm <- gbm %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = VariableFeatures(gbm), npcs = 35, verbose = FALSE)

ElbowPlot(gbm, ndims = 35, reduction = 'pca')

# batch correction on samples
gbm = RunHarmony(gbm, c("orig.ident"), max.iter.harmony = 50, plot_convergence = TRUE)
# harmony converged after 15 iterations

# Inspection of the modalities
p1 = DimPlot_scCustom(gbm, reduction = "harmony", pt.size = 0.1, group.by = "orig.ident")
p2 = DimPlot_scCustom(gbm, reduction = "pca", pt.size = 0.1, group.by = "orig.ident")
pdf(file = 'dimplot_qc_harmony.pdf', width = 14, height = 4)
plot_grid(p2, p1)
dev.off()

p2 = VlnPlot_scCustom(gbm, features = "harmony_1", group.by = "orig.ident",  pt.size = 0)
p1 = VlnPlot_scCustom(gbm, features = "PC_1", group.by = "orig.ident",  pt.size = 0)
pdf(file = 'vlnplot_qc_harmony.pdf', width = 12, height = 4)
plot_grid(p1,p2)
dev.off()

p2 = VlnPlot_scCustom(gbm, features = "harmony_2", group.by = "orig.ident",  pt.size = 0)
p1 = VlnPlot_scCustom(gbm, features = "PC_2", group.by = "orig.ident",  pt.size = 0)
pdf(file = 'vlnplot_qc_harmony2.pdf', width = 12, height = 4)
plot_grid(p1,p2)
dev.off()

# UMAP
gbm <- gbm %>%
  FindNeighbors(reduction = "harmony", dims = 1:35) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(reduction = "harmony", dims = 1:35)

gc()

save_figure('umap', 'seurat_clusters')
save_figure('umap', 'orig.ident')

pdf(file = 'vlnPlot_harmony_comp.pdf', width = 8, height = 6)
VlnPlot(gbm, features = c('harmony_1', 'harmony_2'), group.by = 'orig.ident', pt.size = 0)
dev.off()
pdf(file = 'vlnPlot_pc_comp.pdf', width = 8, height = 6)
VlnPlot(gbm, features = c('PC_1', 'PC_2'), group.by = 'orig.ident', pt.size = 0)
dev.off()

pdf(file = 'featurePlot_qc_1.pdf', width = 9, height = 6)
FeaturePlot(gbm, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'percent.ribo'))
dev.off()
gc()

setwd(seurat_directory)
saveRDS(gbm, 'seu_obj_2310/gbm_harmony.rds')

setwd(seurat_directory)
gbm = readRDS('gbm_harmony.rds')
gbm = FindClusters(gbm, resolution = 0.1)
setwd(fig_directory)
pdf(file = 'umap_res_0.1.pdf', width = 8, height = 6)
DimPlot_scCustom(gbm, reduction = 'umap', group.by = 'seurat_clusters', label = T)
dev.off()

# Cell Annotation -------------------------------------------------------------
## Find Marker Genes ----------------------------------------------------------
# find markers for every cluster compared to all remaining cells, 
# report only the positive ones
gbm = SetIdent(gbm, value = 'RNA_snn_res.0.1')
markers <- FindAllMarkers(gbm, only.pos = TRUE)

setwd(seurat_directory)
saveRDS(markers, file = 'markers_res_01.rds')

markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

pdf(file = 'heatmap_top20_2.pdf', width = 9, height = 20)
DoHeatmap(gbm, features = top20$gene) + 
  scale_fill_gradientn(colours = rev(brewer.pal(11, name = 'PiYG'))) 
dev.off()

saveRDS(gbm, 'gbm_res01.rds')

## Explore Marker Gene Expression from Other Studies --------------------------
# GBmap supplementary tables
# https://www.biorxiv.org/content/10.1101/2022.08.27.505439v1.supplementary-material
gbmap_supp_tab_4 = getGmt('seu_obj_2310/gbmap_supp_table_4.gmt')

# list marker genes of different resolutions of cell annotations based on GBmap
vec_genes = c()
markers_supp_4 = list()
for (i in 1:length(gbmap_supp_tab_4)){
  vec_genes = c(vec_genes, unlist(unname(geneIds(gbmap_supp_tab_4[i]))))
  markers_supp_4 = append(markers_supp_4, 
                          list(unlist(unname(geneIds(gbmap_supp_tab_4[i])))))
}
names(markers_supp_4) = names(gbmap_supp_tab_4)

# Calculate a score for each level of cell annotations
gbm = AddModuleScore(gbm, features = markers_supp_4, name = 'gbmap_cell_anno')

colnames(gbm@meta.data)[40:66] = names(markers_supp_4)

# Visualize score of marker genes of various cell annotation levels
pdf(file = 'featurePlot_gbmap_level_1.pdf', width = 9, height = 3)
FeaturePlot(gbm, features = c('level_1_neoplastic', 'level_1_nonneoplastic'),
            ncol = 2)
dev.off()

pdf(file = 'featurePlot_gbmap_level_2.pdf', width = 9, height = 9)
FeaturePlot(gbm, features = c('level_2_differentiated_like', 'level_2_stem_like',
                              'level_2_glial_neuronal', 'level_2_myeloid',
                              'level_2_lymphoid', 'level_2_vascular'),
            ncol = 2)
dev.off()

pdf(file = 'featurePlot_gbmap_level_3_Neftel.pdf', width = 9, height = 6)
FeaturePlot(gbm, features = c('level_3_AC_like', 'level_3_MES_like',
                              'level_3_OPC_like', 'level_3_NPC_like'),
            ncol = 2)
dev.off()

pdf(file = 'featurePlot_gbmap_level_3_glial_neuronal.pdf', width = 9, height = 9)
FeaturePlot(gbm, features = c('level_3_astrocyte', 'level_3_oligodendrocyte',
                              'level_3_OPC', 'level_3_neuron', 'level_3_RG', 
                              'level_3_endothelial'),
            ncol = 2)
dev.off()

pdf(file = 'featurePlot_gbmap_level_3_myeloid.pdf', width = 9, height = 6)
FeaturePlot(gbm, features = c('level_3_monocyte', 'level_3_TAM_BDM',
                              'level_3_TAM_MG', 'level_3_mast'),
            ncol = 2)
dev.off()

pdf(file = 'featurePlot_gbmap_level_3_lymphoid.pdf', width = 9, height = 9)
FeaturePlot(gbm, features = c('level_3_DC', 'level_3_CD4_CD8',
                              'level_3_NK', 'level_3_B_cell', 'level_3_plasma_B'),
            ncol = 2)
dev.off()


# supplementary table of GBMap for different cell types
# https://www.biorxiv.org/content/10.1101/2022.08.27.505439v1.supplementary-material
gbmap_supp_tables = getGmt('seu_obj_2310/cell_signatures.gmt')

# list marker genes of different cell types based on GBmap
vec_genes = c()
nam_genes = list()
for (i in 1:length(gbmap_supp_tables)){
  vec_genes = c(vec_genes, unlist(unname(geneIds(gbmap_supp_tables[i]))))
  nam_genes = append(nam_genes, list(unlist(unname(geneIds(gbmap_supp_tables[i])))))
}
names(nam_genes) = names(gbmap_supp_tables)
# Calculate a score of marker genes for each cell type
gbm = AddModuleScore(gbm, features = nam_genes, name = 'gbmap_cell_anno')

colnames(gbm@meta.data)[21:39] = names(nam_genes)

# Visualize score of marker genes for various cell types
pdf(file = 'featurePlot_gbmap_malignant.pdf', width = 9, height = 12)
FeaturePlot(gbm, features = c('MES2_like', 'MES1_like', 'AC_like', 'OPC_like',
                              'NPC1_like', 'NPC2_like', 'G1_S', 'G2_M'), 
            ncol = 2)
dev.off()

pdf(file = 'featurePlot_gbmap_brain_resident.pdf', width = 9, height = 6)
FeaturePlot(gbm, features = c('Astrocyte', 'Neuron', 'Oligodendrocyte', 'OPC')
            , ncol = 2)
dev.off()

pdf(file = 'featurePlot_gbmap_vascular.pdf', width = 9, height = 6)
FeaturePlot(gbm, features = c('Endothelial', 'Perivascular_fibroblast', 
                              'Pericyte', 'SMC')
            , ncol = 2)
dev.off()

pdf(file = 'featurePlot_gbmap_immune.pdf', width = 9, height = 6)
FeaturePlot(gbm, features = c('T_cell', 'TAM_MG', 'TAM_BDM')
            , ncol = 2)
dev.off()

gbm = FindClusters(gbm, resolution = 0.5)

pdf(file = 'dimplot_res05.pdf', width = 8, height = 6)
DimPlot_scCustom(gbm, reduction = 'umap')
dev.off()

## Identify malignant cells ---------------------------------------------------
mal = subset(gbm, idents = c('1', '2', '3', '4', '5', '6', '7', '8', '10', '11', 
                             '12', '13', '20', '21', '22', '24', '29', '33'))

### Neftel's glioma cell states -----------------------------------------------
# https://www.sciencedirect.com/science/article/pii/S0092867419306877
neftel_supp_tab_2 = getGmt('seu_obj_2310/malignant_signature_Neftel.gmt')

# list marker genes of Neftel's glioma cell states
vec_genes = c()
nam_genes = list()
for (i in 1:length(neftel_supp_tab_2)){
  vec_genes = c(vec_genes, unlist(unname(geneIds(neftel_supp_tab_2[i]))))
  nam_genes = append(nam_genes, list(unlist(unname(unique(geneIds(neftel_supp_tab_2[i]))))))
}
names(nam_genes) = names(neftel_supp_tab_2)

# Calculate a score of marker genes for each glioma state
mal = AddModuleScore(mal, features = nam_genes, name = 'Neftel')
colnames(mal@meta.data)[68:72] = names(nam_genes)

# Visualize score of Neftel's marker genes for each glioma subtype
pdf(file = 'mal_neftel_featureplot.pdf', width = 9, height = 9)
FeaturePlot(mal, features = c('MES_like_Neftel', 'AC_like_Neftel', 
                              'OPC_like_Neftel', 'NPC_like_Neftel', 
                              'Cell_Cycle_Neftel'), ncol = 2)
dev.off()

# Prepare scores for butterfly plot
df = as.data.frame(matrix(ncol = 4, nrow = dim(mal)[2]))
row.names(df) = mal@assays$RNA@data@Dimnames[[2]]
df$V1 = mal$AC_like_Neftel
df$V2 = mal$OPC_like_Neftel
df$V3 = mal$MES_like_Neftel
df$V4 = mal$NPC_like_Neftel

df$V1 = as.double(df$V1)
df$V2 = as.double(df$V2)
df$V3 = as.double(df$V3)
df$V4 = as.double(df$V4)

ho = hierarchy(df) # from scrabble package
gr = mal$RNA_snn_res.0.1
gr_val = unname(gr)
gr = names(gr)
names(gr) = gr_val
data("palette36")

pal = palette36[as.numeric(levels(mal$RNA_snn_res.0.1))+1]
names(pal) = levels(mal$RNA_snn_res.0.1)
pal  

# Visualize butterfly plot
pdf(file = 'butterflyplot_Neftel_res01.pdf', width = 9, height = 8)
plot_hierarchy(ho, quadrant.names = c('AC-like', 'OPC-like', 'MES-like', 'NPC-like'), 
               groups = gr, group.cols = pal)
dev.off()

# Assign Neftel's glioma states based on the quadrant a cell is placed in a butterfly plot
df2 = as.data.frame(matrix(ncol = 1, nrow = dim(mal)[2]))
row.names(df2) = mal@assays$RNA@data@Dimnames[[2]]
for (i in row.names(ho)){
  if (ho[i, 'X'] <= 0 & ho[i, 'Y'] <= 0){
    df2[i, 'V1'] = 'AC-like'
  } else if (ho[i, 'X'] <= 0 & ho[i, 'Y'] >= 0){
    df2[i, 'V1'] = 'MES-like'
  } else if (ho[i, 'X'] >= 0 & ho[i, 'Y'] >= 0){
    df2[i, 'V1'] = 'NPC-like'
  } else if (ho[i, 'X'] >= 0 & ho[i, 'Y'] <= 0){
    df2[i, 'V1'] = 'OPC-like'
  }
}

# Add Neftel's cell states to seurat object
mal = AddMetaData(mal, df2$V1, col.name = 'malignant_Neftel')

## Assign clusters to cell types ----------------------------------------------
gbm@meta.data[['my_cell_anno']] = NA
df3 = data.frame(gbm@meta.data[['my_cell_anno']], row.names = rownames(gbm@meta.data))
colnames(df3) = c('my_cell_anno')
df3[rownames(df2), 'my_cell_anno'] = df2$V1

res05 = data.frame(gbm$RNA_snn_res.0.5)
my_cell_anno = cbind(res05, df3)

my_cell_anno[my_cell_anno$gbm.RNA_snn_res.0.5 == 25, 'my_cell_anno'] = 'Endothelial'
my_cell_anno[my_cell_anno$gbm.RNA_snn_res.0.5 == 28, 'my_cell_anno'] = 'Pericyte'
my_cell_anno[(my_cell_anno$gbm.RNA_snn_res.0.5 == 27) |
               (my_cell_anno$gbm.RNA_snn_res.0.5 == 26) |
               (my_cell_anno$gbm.RNA_snn_res.0.5 == 9) |
               (my_cell_anno$gbm.RNA_snn_res.0.5 == 0) |
               (my_cell_anno$gbm.RNA_snn_res.0.5 == 17) |
               (my_cell_anno$gbm.RNA_snn_res.0.5 == 30) |
               (my_cell_anno$gbm.RNA_snn_res.0.5 == 19) |
               (my_cell_anno$gbm.RNA_snn_res.0.5 == 15) |
               (my_cell_anno$gbm.RNA_snn_res.0.5 == 18)
             , 'my_cell_anno'] = 'Myeloid'
my_cell_anno[my_cell_anno$gbm.RNA_snn_res.0.5 == 14, 'my_cell_anno'] = 'Lymphoid'
my_cell_anno[(my_cell_anno$gbm.RNA_snn_res.0.5 == 16) |
               (my_cell_anno$gbm.RNA_snn_res.0.5 == 23) |
               (my_cell_anno$gbm.RNA_snn_res.0.5 == 31) |
               (my_cell_anno$gbm.RNA_snn_res.0.5 == 32) 
             , 'my_cell_anno'] = 'Oligo'
any(is.na(my_cell_anno$my_cell_anno))

# save cell types in metadata
gbm@meta.data[['my_cell_anno']] = my_cell_anno$my_cell_anno

gbm = SetIdent(gbm, value = 'my_cell_anno')

pdf(file = 'dimplot_gbm_my_cell_anno2.pdf', width = 8, height = 6)
DimPlot_scCustom(gbm, reduction = 'umap', pt.size = 0.01, label = F)
dev.off()

saveRDS(gbm, 'seu_obj_2310/gbm_my_cell_anno.rds')

# Expression of marker genes in defined cell types
pdf(file = 'dotplot_gbm_combined_markers.pdf', width = 6, height = 7)
Clustered_DotPlot(gbm, features = c('VWF', 'MZB1', 'MS4A1', 'NKG7', 'CD3D', 'CX3CR1',
                                    'TGFBI', 'GAD1', 'PDGFRB', 'ACTA2', 'PECAM1', 
                                    'P2RY12', 'GFAP', 'PDGFRA', 'EGFR', 'VIM', 'PTPRC',
                                    'MT1X', 'BCAN', 'CDK4', 'AQP4', 'VCAN', 'TEK',
                                    'MOBP', 'MOG', 'PLP1', 'MBP'),
                  row_label_size = 12, plot_km_elbow = F, show_parent_dend_line = F,
                  colors_use_exp = colorRampPalette(rev(brewer.pal(11, name = 'RdYlBu')))(12)) 

dev.off()

# Proportion plots
plot_integrated_clusters(gbm, 'seurat_clusters', 'my_cell_anno', 'res05_cell_anno')
plot_integrated_clusters(gbm, 'seurat_clusters', 'orig.ident', 'res05_patient_sample')
plot_integrated_clusters(gbm, 'orig.ident', 'my_cell_anno', is_cluster_x_int = F, 
                         'patient_sample_cell_anno')
plot_integrated_clusters(gbm, 'my_cell_anno', 'orig.ident', is_cluster_x_int = F, 
                         'cell_anno_patient_sample')

plot_integrated_clusters(mal, 'RNA_snn_res.0.1', 'malignant_Neftel', 
                         'mal_res01_cell_anno')

# Expression of Neftel's glioma markers in glioma subtypes
mal = SetIdent(mal, value = 'malignant_Neftel')
pdf(file = 'dotplot_mal_neftel_cellcycle_markers.pdf', width = 5, height = 11)
Clustered_DotPlot(mal, features = unique(nam_genes$Cell_Cycle_Neftel),
                  row_label_size = 12, plot_km_elbow = F, show_parent_dend_line = F,
                  colors_use_exp = colorRampPalette(rev(brewer.pal(11, name = 'RdYlBu')))(12)) 
dev.off()

# Expression of proliferative markers (cell cycle state) in clusters 
Clustered_DotPlot(gbm, features = c('TUBB', 'USP1', 'EZH2', 'TOP2A', 
                                    'MKI67', 'AURKB', 'FOXM1', 'TYMS'), 
                  row_label_size = 12, plot_km_elbow = F, show_parent_dend_line = F,
                  colors_use_exp = colorRampPalette(rev(brewer.pal(11, name = 'RdYlBu')))(12)) 


my_metadata = data.frame(matrix(nrow = ncol(gbm), ncol = 5))
rownames(my_metadata) = colnames(gbm)
colnames(my_metadata) = c('sample', 'patient', 'sampleSite', 'site', 'sampleSiteCellAnno')
my_metadata$sample = gbm$orig.ident

# Annotate cell cycle cluster
my_metadata[['rna_res_05']] = gbm$RNA_snn_res.0.5
my_metadata[my_metadata$rna_res_05 == '7', 'cell_anno_2'] = 'CellCycle'

# Add clinical metadata to seurat object
samples = strsplit(gbm$orig.ident, split = '_')
samples = lapply(samples, function(x) paste(x[2], x[3], sep = ''))
donor_id = unlist(unname(samples))

my_metadata$patient = donor_id

my_metadata$site = c(rep('NonEnhancing', ncol(gbm)))

my_metadata[(my_metadata$sample == 'Sample_4554_DH_2') |
              (my_metadata$sample == 'Sample_4554_DH_3') |
              (my_metadata$sample == 'Sample_5363_SR_1') |
              (my_metadata$sample == 'Sample_5363_SR_2') |
              (my_metadata$sample == 'Sample_5675_SS_1') |
              (my_metadata$sample == 'Sample_6025_SS_1') |
              (my_metadata$sample == 'Sample_6369_SS_2'), 'site'] = 'Enhancing'
my_metadata$sampleSite = paste(my_metadata$sample, my_metadata$site, sep = '')

my_metadata[['cell_anno_2']] = gbm$my_cell_anno
my_metadata$sampleSiteCellAnno = paste(my_metadata$sampleSite, gbm$my_cell_anno, sep = '')

gbm = AddMetaData(gbm, my_metadata)

gbm = SetIdent(gbm, value = 'cell_anno_2')

# Proportion plots
plot_integrated_clusters(gbm, 'cell_anno_2', 'sampleSite', F, 'cell_anno_2_sample')

# Pie plot - cell type proportions in the whole dataset
prop = table(as.character(gbm$cell_anno_2))
no_all_cells = sum(prop)
prop = prop / no_all_cells * 100
prop = as.data.frame(prop)
rownames(prop) = prop$Var1
ordered_name = c('AC-like', 'NPC-like', 'OPC-like', 'MES-like', 'CellCycle',
                 'Myeloid', 'Oligo', 'Endothelial', 'Pericyte', 'Lymphoid')
prop = prop[ordered_name, ]

prop
pdf(file = 'pie_cell_anno_2.pdf', height = 3, width = 3)
pie(prop$Freq , labels = prop$Var1, border = "black", col = my_colors)
dev.off()

prop['AC-like', 'Freq'] + prop['NPC-like', 'Freq'] + prop['OPC-like', 'Freq'] + prop['MES-like', 'Freq'] + prop['CellCycle', 'Freq']

gbm$cell_anno_2 = factor(gbm$cell_anno_2, levels = c('AC-like', 'NPC-like', 
                                                     'OPC-like', 'MES-like', 
                                                     'CellCycle',
                                                     'Myeloid', 'Oligo',
                                                     'Endothelial', 'Pericyte', 
                                                     'Lymphoid'))

saveRDS(gbm, file = 'seu_obj_2310/gbm_w_metadata.rds')
