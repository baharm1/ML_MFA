# Import libraries ------------------------------------------------------------
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(ggbeeswarm)
library(ggpubr)

# Plot CNN-predicted relative de novo GMP synthesis flux MMF-treated mice ----

mice_pred_dir = './gmp_denovo_glioma_config_20240517_1138/mice_mmf_pred_gmp_denovo_glioma'
mice_pred_files = list.files(path = mice_pred_dir)

list_mice_pred = c()
for (ps in 1:length(mice_pred_files)){
  mice_pred = read.delim(file = paste(mice_pred_dir, '/', mice_pred_files[ps], 
                                      sep = ''), 
                         sep = '\t', header = T)
  colnames(mice_pred) = c('ind', 'denovo')
  list_mice_pred = c(list_mice_pred, list(mice_pred$denovo))
}

names(list_mice_pred) = substr(mice_pred_files, 1, 10)

mice_pred_df = dplyr::bind_rows(list_mice_pred)

drops <- c("high2_773", "high2_774", "high2_780",
           "high4_775", "high4_778", "high4_779",
           "low2_766G", "low2_767G", "low2_782G",
           "low4_769G", "low4_771G", "low4_772G",
           "high1_764", "high1_765", "high1_781", )

drops <- c("high2_773G", "high2_774G", "high2_780G",
           "high4_775G", "high4_778G", "high4_779G")

mice_pred_df = mice_pred_df[, !(names(mice_pred_df) %in% drops)]


mice_pred_long = tidyr::pivot_longer(mice_pred_df, cols = 1:ncol(mice_pred_df))

treatment_sample = stringr::str_split_fixed(mice_pred_long$name, "_", 2)

treatment_sample = as.data.frame(treatment_sample)
mice_pred_long[['treatment']] = treatment_sample$V1
mice_pred_long[['sample']] = treatment_sample$V2

ReplicateAverages = mice_pred_long %>% 
  group_by(treatment, sample) %>% 
  summarise_each(list(mean)) 

TreatmentAverages = mice_pred_long %>% 
  group_by(treatment) %>% 
  summarise_each(list(mean)) 

ReplicateAverages.summary <- ReplicateAverages %>%
  group_by(treatment) %>%
  summarise(sd = sd(value, na.rm = T),
            value = mean(value))

ReplicateAverages.summary

p = ggplot(ReplicateAverages, aes(treatment, value)) +
  geom_bar(stat = "identity", data = ReplicateAverages.summary,
           fill = NA, color = 'black') +
  geom_jitter(position = position_jitter(0.2),
              color = 'black') +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),
                data = ReplicateAverages.summary, width = 0.2)

p
#https://rupress.org/jcb/article/219/6/e202001064/151717/SuperPlots-Communicating-reproducibility-and
mice_pred_long_2 = mice_pred_long[, -c(1)]
mice_pred_long_3 = mice_pred_long_2[, -c(3)]

mice_pred_long_3$treatment = factor(mice_pred_long_3$treatment, 
                                    levels = c('ctrl', 'low1', 'high1', 
                                               'low2', 'low4'))

mice_pred_long_2$treatment = factor(mice_pred_long_2$treatment, 
                                    levels = c('ctrl', 'low1', 'high1',
                                               'low2', 'low4'))

ReplicateAverages = mice_pred_long_2 %>% 
  group_by(treatment, sample) %>% 
  summarise_each(list(mean)) 

TreatmentAverages = mice_pred_long_3 %>% 
  group_by(treatment) %>% 
  summarise_each(list(mean)) 

#ReplicateAverages[['sample']] = NULL

p = ggplot(mice_pred_long_2, aes(x=treatment,y=value,color=factor(sample))) + 
  #geom_jitter() +
  geom_beeswarm(corral = "wrap", corral.width = 0.9, cex = 1, size = 0.5) + 
  #scale_colour_brewer(palette = "Set3") + 
  scale_color_manual(values = c('#90AD1C', #CTRL
                                '#C075A6', #LOW1
                                '#AA0DFE', #LOW1
                                '#1C8356', #CTRL
                                '#DEA0FD', #LOW1
                                '#325A9B', #HIGH1
                                '#2ED9FF', #HIGH1
                                '#C4451C', #LOW2
                                '#FEAF16', #LOW2 
                                '#822E1C', #LOW4
                                '#3283FE', #LOW4
                                '#3a002c', #LOW4
                                '#600049', #HIGH1
                                '#b080a4'  #LOW2
  )) +
  stat_compare_means(data = mice_pred_long_2,
                     comparisons = list(c('ctrl', 'low1'),
                                        c('ctrl', 'high1'),
                                        c('ctrl', 'low2'),
                                        c('ctrl', 'low4'),
                                        c('low1', 'high1'),
                                        c('low1', 'low2'),
                                        c('low1', 'low4'),
                                        c('low2', 'low4')),
                     #method = 'anova',
                     label.y = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.95),
                     size = 2) + 
  geom_beeswarm(data = ReplicateAverages, size = 3, shape = 15) + 
  geom_beeswarm(data = TreatmentAverages, size = 3, shape = 15, color = '#474747') +
  theme_classic()
#theme(legend.position="none")

ggsave(plot = p, filename = 'mice_mmf_glioma_240517_1138_gmp_3.pdf',
       device = "pdf", height = 5, #5.5
       width = 6, units = "in")#3.2


p = ggplot(mice_pred_long_2, aes(x=treatment,y=value,color=factor(sample))) + 
  #geom_jitter() +
  geom_beeswarm(corral = "wrap", corral.width = 0.9, cex = 1, size = 0.5) + 
  #scale_colour_brewer(palette = "Set3") + 
  scale_color_manual(values = c('#90AD1C', #CTRL
                                '#C075A6', #LOW1
                                '#AA0DFE', #LOW1
                                '#1C8356', #CTRL
                                '#DEA0FD', #LOW1
                                '#C4451C', #high1
                                '#FEAF16', #high1
                                '#822E1C' #high1
  )) +
  stat_compare_means(#data = ReplicateAverages,
    comparisons = list(c('ctrl', 'low1'),
                       c('ctrl', 'high1'),
                       c('low1', 'high1')),
    #method = 'anova',
    label.y = c(0.05, 0.1, 0.95),
    size = 2) + 
  geom_beeswarm(data = ReplicateAverages, size = 3, shape = 15) + 
  geom_beeswarm(data = TreatmentAverages, size = 3, shape = 15, color = '#474747') +
  theme_classic()
#theme(legend.position="none")

ggsave(plot = p, filename = 'mice_mmf_glioma_240517_1138_imputed_5e4_rm0_100_conc.pdf',
       device = "pdf", height = 5, #5.5
       width = 6, units = "in")#3.2



# 2406
names(list_mice_pred) = c('ctrl', 'low1', 'low2', 'low4')

mice_pred_df = dplyr::bind_rows(list_mice_pred)

mice_pred_long = tidyr::pivot_longer(mice_pred_df, cols = 1:ncol(mice_pred_df))

df.summary <- mice_pred_long %>%
  group_by(name) %>%
  summarise(sd = sd(value, na.rm = T),
            value = mean(value))
df.summary

p = ggplot(mice_pred_long, aes(name, value)) +
  geom_bar(stat = "identity", data = df.summary,
           fill = NA, color = 'black') +
  geom_jitter(position = position_jitter(0.2),
              color = 'black') +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),
                data = df.summary, width = 0.2)

p

p = ggplot(mice_pred_long, aes(name, value)) +
  geom_violin(trim = F, scale = 'width')


a1 = aov(data = ReplicateAverages, formula = value ~ treatment)
pairwise.t.test(ReplicateAverages$value, ReplicateAverages$treatment, 
                p.adj = "bonf")
summary(a1)
TukeyHSD(a1)

# GBM12, GBM38, HF2303 --------------------------------------------------------
mouse_pred_dir  = './serine_plasma_glioma_config_20240408_2102/mice_pred_serine_plasma_glioma/'
pred_nat_rev = list.files(path = mouse_pred_dir, recursive = T)
mouse_pred = list()

for (i in pred_nat_rev){
  mouse_file = read.delim(file = paste(mouse_pred_dir, i, sep = ''), 
                          header = T, quote = "")
  mouse_cond = substr(i, 17, nchar(i) - 11)
  mouse_pred[[mouse_cond]] = mouse_file$pred_serine_plasma_glioma_config_mice_nat_rev
  
}
mouse_pred = as.data.frame(mouse_pred)

mouse_df = data.frame('sample' = c(rep('GBM12', 200), rep('GBM38', 200), 
                                       rep('HF2303', 200)),
                      'cond' = c(rep('Tumor', 100), rep('Tumor-SG', 100),
                                 rep('Tumor', 100), rep('Tumor-SG', 100),
                                 rep('Tumor', 100), rep('Tumor-SG', 100)),
                      'value' = c(mouse_pred$GBM12, mouse_pred$GBM12_SG,
                                  mouse_pred$GBM38, mouse_pred$GBM38_SG,
                                  mouse_pred$HF2303, mouse_pred$HF2303_SG))

                                       

p = ggplot(mouse_df, aes(x = sample, y = value, fill = cond)) +
  geom_violin(trim = T, scale = 'width', position = position_dodge(0.8), 
              alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.8)) +
  stat_summary(fun=mean, geom="point", size=1, shape = 1, 
               position = position_dodge(0.8), color="black") + 
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_fill_manual(values=c("#f16623", "#23aef1")) 

p

ggsave(filename = 'mice_pred_serine_plasma_glioma_config_mice_nat_rev.pdf',
       plot = p, device = 'pdf', width = 4.5, height = 3)

stats_plasma_ser = ggpubr::compare_means(data = mouse_df, 
                                       formula = value ~ cond, 
                                       group.by = 'sample', 
                                       p.adjust.method = 'holm')
write.csv(stats_plasma_ser, file = "stats_mice_pred_serine_plasma_glioma_config_mice_nat_rev_cond_wilcox_holm.csv")

stats_plasma_ser = ggpubr::compare_means(data = mouse_df, 
                                         formula = value ~ sample, 
                                         group.by = 'cond', 
                                         p.adjust.method = 'holm')

write.csv(stats_plasma_ser, file = "stats_mice_pred_serine_plasma_glioma_config_mice_nat_rev_wilcox_holm.csv")

# TRP vs. GBM38 ---------------------------------------------------------------
mouse_pred_dir  = './gmp_denovo_glioma_config_20240517_1138/mice_pred_gmp_denovo_glioma_TRP_GBM38/'
pred_nat_rev = list.files(path = mouse_pred_dir, recursive = T)
mouse_pred = list()

for (i in pred_nat_rev){
  mouse_file = read.delim(file = paste(mouse_pred_dir, i, sep = ''), 
                          header = T, quote = "")
  mouse_cond = substr(i, 1, nchar(i) - 23)
  mouse_pred[[mouse_cond]] = mouse_file$pred_gmp_denovo_TRP_GBM38_imputed_config
  
}
mouse_pred = as.data.frame(mouse_pred)

mouse_df = data.frame('sample' = c(rep('TRP', 100), rep('GBM38', 100)),
                      'value' = c(mouse_pred$TRP, mouse_pred$GBM38))
mouse_df$sample = factor(mouse_df$sample, levels = c('TRP', 'GBM38'))

p = ggplot(mouse_df, aes(x = sample, y = value, fill = sample)) +
  geom_violin(trim = T, scale = 'width', position = position_dodge(0.8), 
              alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.8)) +
  stat_summary(fun=mean, geom="point", size=1, shape = 1, 
               position = position_dodge(0.8), color="black") + 
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_fill_manual(values=c("#7ecbec", "#f795bd")) 

p

ggsave(filename = 'mice_pred_gmp_denovo_TRP_GBM38.pdf',
       plot = p, device = 'pdf', width = 2, height = 3)

stats_gmp = ggpubr::compare_means(data = mouse_df, 
                                  formula = value ~ sample, 
                                  p.adjust.method = 'holm')

write.csv(stats_gmp, file = "stats_mice_pred_gmp_denovo_TRP_GBM38_holm.csv")


# correlation plots -----------------------------------------------------------
# Figure 3I, S3L-M 
# plot M+3 serine to M+3 PG and M+1 serine to M+1 plasma serine
ML_pred_dir = "./all_patients_pred_serine_CNN/"

glioma_denovo = read.delim(file = paste(ML_pred_dir, 'pred_denovo_glioma_20240408_1327.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')

glioma_plasma = read.delim(file = paste(ML_pred_dir, 'pred_plasma_glioma_20240408_2102.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
glioma_tme = read.delim(file = paste(ML_pred_dir, 'pred_tme_glioma_20240408_1327_2102.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
cortex_denovo = read.delim(file = paste(ML_pred_dir, 'pred_denovo_cortex_20240508_1220.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
cortex_plasma = read.delim(file = paste(ML_pred_dir, 'pred_plasma_cortex_20240508_1220.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')

patient_site_names = colnames(glioma_denovo)

patient_folders = list.dirs(path = "./data/patient_data/patient_mid_mc_serine/",
                            recursive = F)
mid_names = read.delim(file = "./data/patient_data/mid_name_patient_serine.txt", 
                       header = F, sep = '\t')

for (ps in 1:ncol(glioma_denovo)){
patient_mids = read.delim(file = patient_folders[ps], header = F, sep = '\t')

colnames(patient_mids) = mid_names$V1

glioma_denovo[[paste(patient_site_names[ps], 'M3', sep = '_')]] = patient_mids$SER3 / patient_mids$PG3
glioma_plasma[[paste(patient_site_names[ps], 'M1', sep = '_')]] = patient_mids$SER1 / patient_mids$SERp1

glioma_tme[[paste(patient_site_names[ps], 'M1', sep = '_')]] = patient_mids$SER1 / patient_mids$SERc1
glioma_tme[[paste(patient_site_names[ps], 'M2', sep = '_')]] = patient_mids$SER2 / patient_mids$SERc2
glioma_tme[[paste(patient_site_names[ps], 'M3', sep = '_')]] = patient_mids$SER3 / patient_mids$SERc3

cortex_denovo[[paste(patient_site_names[ps], 'M3', sep = '_')]] = patient_mids$SERc3 / patient_mids$PGc3
cortex_plasma[[paste(patient_site_names[ps], 'M1', sep = '_')]] = patient_mids$SERc1 / patient_mids$SERp1

}

# Figure S3K-L
glioma_denovo_cor = calculateCorrelation(glioma_denovo)
pdf(file = 'corrplot_glioma_denovo.pdf', width = 10, height = 10)
p = correlationPlot(glioma_denovo_cor$M, glioma_denovo_cor$adj_pval)
print(p)
dev.off()

glioma_plasma_cor = calculateCorrelation(glioma_plasma)
pdf(file = 'corrplot_glioma_plasma.pdf', width = 10, height = 10)
p = correlationPlot(glioma_plasma_cor$M, glioma_plasma_cor$adj_pval)
print(p)
dev.off()

glioma_tme_cor = calculateCorrelation(glioma_tme)
pdf(file = 'corrplot_glioma_tme.pdf', width = 20, height = 20)
p = correlationPlot(glioma_tme_cor$M, glioma_tme_cor$adj_pval)
print(p)
dev.off()

glioma_denovo_ps = colMeans(glioma_denovo)
glioma_plasma_ps = colMeans(glioma_plasma)

cortex_denovo_ps = colMeans(cortex_denovo)
cortex_plasma_ps = colMeans(cortex_plasma)

glioma_ps = rbind(glioma_denovo_ps[16:30], 
                  glioma_plasma_ps[16:30],
                  cortex_denovo_ps[16:30],
                  cortex_plasma_ps[16:30])

glioma_ps_ML = rbind(glioma_denovo_ps[1:15],
                    glioma_plasma_ps[1:15],
                    cortex_denovo_ps[1:15],
                    cortex_plasma_ps[1:15])

cortex_ps_ML = rbind(cortex_denovo_ps,
                     cortex_plasma_ps)

rownames(cortex_ps_ML) = c('denovo_c', 'plasma_c')

rownames(glioma_ps) = c('SERg3/PGg3', 
                        'SERg1/SERp1', 
                        'SERc3/PGc3', 
                        'SERc1/SERp1')

rownames(glioma_ps_ML) = c('denovo_g', 
                           'plasma_g', 
                           'denovo_c', 
                           'plasma_c')

scaled_glioma_ps = scale(t(glioma_ps))
scaled_glioma_ps_ML = scale(t(glioma_ps_ML))
scaled_cortex_ps_ML = scale(t(cortex_ps_ML))
  
pdf(file = 'pheatmap_patient_MID_ratio.pdf', width = 10, height = 3)
pheatmap(t(scaled_glioma_ps), cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = COL1('OrRd'))
dev.off()

pdf(file = 'pheatmap_patient_ML_ratio_OrRd.pdf', width = 10, height = 3)
pheatmap(t(scaled_glioma_ps_ML), cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = COL1('OrRd'))#YlGn
dev.off()

cortex_denovo = read.delim(file = paste(ML_pred_dir, 'ML_denovo_cortex_1615_NandE.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
cortex_plasma = read.delim(file = paste(ML_pred_dir, 'ML_plasma_cortex_1615_NandE.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')

pdf(file = 'pheatmap_patient_ML_ratio_cortex_OrRd.pdf', width = 5, height = 3)
pheatmap(t(scaled_cortex_ps_ML), cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = COL1('OrRd'))
dev.off()

plasma_all = read.delim(file = paste(ML_pred_dir, 'ML_plasma_1327_2102_1615_NandE.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
scaled_denovo = scale(t(colMeans(denovo_all)))
pheatmap(t(scaled_denovo), cluster_rows = F, cluster_cols = F, show_rownames = T,
         color = COL1('YlGn'))

# Figure 3I -------------------------------------------------------------------
#240523 - scatter plots - correlation of scaled flux & MID ratio
scaled_glioma_ps = data.frame(scaled_glioma_ps)
scaled_glioma_ps_ML = data.frame(scaled_glioma_ps_ML)

denovog_ser3 = data.frame('SERg3.PG3' = scaled_glioma_ps$SERg3.PGg3,
                          'denovo_g' = scaled_glioma_ps_ML$denovo_g,
                          'patient' = factor(seq(1, nrow(scaled_glioma_ps))))

plasmag_ser1 = data.frame('SERg1.SERp1' = scaled_glioma_ps$SERg1.SERp1,
                          'plasma_g' = scaled_glioma_ps_ML$plasma_g,
                          'patient' = factor(seq(1, nrow(scaled_glioma_ps))))
library("ggpubr")
p = ggscatter(denovog_ser3, x = "SERg3.PG3", y = "denovo_g", color = "patient",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "SERg3/PGg3", ylab = "Glucose-derived serine synthesis flux",
          label = "patient") +
  geom_smooth(method = "lm", color = "black") 
  
ggsave(filename = paste('scatter', 'denovo_c_SERc3_PGc3_cor',
                        '.pdf', sep = '_'), 
       plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')

p = ggscatter(plasmag_ser1, x = "SERg1.SERp1", y = "plasma_g", color = "patient",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "SERg1/SERp1", ylab = "Plasma serine uptake flux",
              label = "patient") +
  geom_smooth(method = "lm", color = "black") 

ggsave(filename = paste('scatter', 'plasma_c_SERc1_SERp1_cor',
                        '.pdf', sep = '_'), 
       plot = p, device = 'pdf', width = 5, height = 5.5, units = 'in')

# add cortex MIDs from nonenhancing and enhancing 
cortex_denovo = read.delim(file = paste(ML_pred_dir, 'ML_denovo_cortex_1615_NandE.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')
cortex_plasma = read.delim(file = paste(ML_pred_dir, 'ML_plasma_cortex_1615_NandE.txt', 
                                        sep = ''), 
                           header = T, sep = '\t')

for (ps in 1:(ncol(cortex_denovo)-1)){
  patient_mids = read.delim(file = paste(patient_folders[2*ps-1], '/mid_mc.txt', sep = ''), 
                            header = F, sep = '\t')
  
  colnames(patient_mids) = mid_names$V1
  E3 = patient_mids$SERc3 / patient_mids$PGc3
  E1 = patient_mids$SERc1 / patient_mids$SERp1
  
  patient_mids = read.delim(file = paste(patient_folders[2*ps], '/mid_mc.txt', sep = ''), 
                            header = F, sep = '\t')
  
  colnames(patient_mids) = mid_names$V1
  N3 = patient_mids$SERc3 / patient_mids$PGc3
  N1 = patient_mids$SERc1 / patient_mids$SERp1
  
  cortex_denovo[[paste(patient_site_names[2*ps-1], 'M3', sep = '_')]] = c(E3, N3)
  cortex_plasma[[paste(patient_site_names[2*ps-1], 'M1', sep = '_')]] = c(E1, N1)
  
}
ps = 8
patient_mids = read.delim(file = paste(patient_folders[2*ps-1], '/mid_mc.txt', sep = ''), 
                          header = F, sep = '\t')

colnames(patient_mids) = mid_names$V1
E3 = patient_mids$SERc3 / patient_mids$PGc3
E1 = patient_mids$SERc1 / patient_mids$SERp1

cortex_denovo[[paste(patient_site_names[2*ps-1], 'M3', sep = '_')]] = c(E3, E3)
cortex_plasma[[paste(patient_site_names[2*ps-1], 'M1', sep = '_')]] = c(E1, E1)

library(RcmdrMisc)
cortex_denovo_cor = calculateCorrelation(cortex_denovo)
pdf(file = 'corrplot_cortex_denovo_200.pdf', width = 10, height = 10)
p = correlationPlot(cortex_denovo_cor$M, cortex_denovo_cor$adj_pval)
print(p)
dev.off()

cortex_plasma_cor = calculateCorrelation(cortex_plasma)
pdf(file = 'corrplot_cortex_plasma_200.pdf', width = 10, height = 10)
p = correlationPlot(cortex_plasma_cor$M, cortex_plasma_cor$adj_pval)
print(p)
dev.off()


cortex_denovo_ps = colMeans(cortex_denovo)
cortex_plasma_ps = colMeans(cortex_plasma)

cortex_ps = rbind(cortex_denovo_ps[9:16],
                  cortex_plasma_ps[9:16])

cortex_ps_ML = rbind(cortex_denovo_ps[1:8],
                     cortex_plasma_ps[1:8])

rownames(cortex_ps) = c('SERg3/PGg3', 
                        'SERg1/SERp1')

rownames(cortex_ps_ML) = c('denovo_g', 
                           'plasma_g')

scaled_cortex_ps = scale(t(cortex_ps))
scaled_cortex_ps_ML = scale(t(cortex_ps_ML))

scaled_glioma_ps = scaled_cortex_ps
scaled_glioma_ps_ML = scaled_cortex_ps_ML

# functions -------------------------------------------------------------------

correlationPlot = function(M, adj_pval){
  p = corrplot(M, type = "full", method = 'ellipse', col = rev(COL2('PiYG')),
               tl.col = "black", tl.cex = 1.5, p.mat = adj_pval, 
               sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.1,
               insig = 'label_sig', diag = F) 
  return(p)
}

calculateCorrelation = function(df){
  M = cor(df, method = 'pearson')
  M = M[rowSums(is.na(M)) != ncol(M) - 1, colSums(is.na(M)) != nrow(M) - 1]
  
  adj_pval = rcorr.adjust(df)
  adj_pval = adj_pval$P
  adj_pval = ifelse(adj_pval == "<.0001", 0.0001, as.numeric(adj_pval))
  adj_pval = adj_pval[rowSums(is.na(adj_pval)) != ncol(adj_pval), 
                      colSums(is.na(adj_pval)) != nrow(adj_pval)]
  
  return(list('M' = M, 'adj_pval' = adj_pval))
}
