# MEBOCOST metabolite connectivity between cells

library(ggplot2)

commu_df = read.csv(file = './darmanis_MEBOCOST_output/commu_df.csv',
                    header = T, quote = '')

commu_df = commu_df[commu_df$permutation_test_fdr < 0.05, ]

metab_connectivity = as.data.frame(table(commu_df$Metabolite_Name))

metab_connectivity = metab_connectivity[order(metab_connectivity$Freq,
                                              decreasing = F),]
metab_connectivity$Var1 = factor(metab_connectivity$Var1, 
                                 levels = metab_connectivity$Var1)

p = ggplot(metab_connectivity, aes(x = Var1, y = Freq, fill = Var1))+
  geom_bar(width = 0.7, stat = "identity")+
  coord_flip() +
  theme_classic() +
  theme(legend.position="none") 

ggsave(filename = 'MEBOCOST_sig_metabolite_connectivity.pdf', plot = p, 
       device = 'pdf', height = 5, width = 4, units = 'in')

