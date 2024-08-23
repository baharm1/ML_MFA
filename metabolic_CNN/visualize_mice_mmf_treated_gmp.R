# Import libraries ------------------------------------------------------------
library(tidyr)
library(RColorBrewer)
library(see)
library(ggplot2)
library(dplyr)
library(gghalves)
library(stringr)
library(Seurat)
library(ggbeeswarm)
library(ggpubr)

# Plot CNN-predicted relative de novo GMP synthesis flux MMF-treated mice ----

mice_pred_dir = './gmp_denovo_glioma_config_20240517_1138/mice_pred_gmp_denovo_glioma_config'
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

