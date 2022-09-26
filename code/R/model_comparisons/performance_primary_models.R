library(tidyverse)

df <- data.frame(Data = factor(c("clinical", "PRS", "transcriptomic", "proteomic", "lipids and \nmetabolites"), 
                               levels = c("clinical", "PRS", "transcriptomic", "proteomic", "lipids and \nmetabolites")), 
                 Accuracy = c(0.62, 0.59, 0.59, 0.67, 0.55),
                 AUROC = c(0.63, 0.61, 0.57, 0.67, 0.56))

plot_unimodal <- ggplot(df, aes(Data, AUROC)) + 
  geom_point(aes(colour=AUROC), size = 5) +
  scale_color_continuous(low = "yellow", high = "red", na.value = NA, limits=c(0.5, 0.8)) +
  theme_light() +
  ylim(0.5, 1.0) +
  geom_text(aes(label = format(AUROC, nsmall = 2), vjust = AUROC-2), color = "black", size = 6) +
  theme(text = element_text(size=24),
        axis.text=element_text(colour="black"),
        legend.justification = 'top') +
  scale_x_discrete(limits=rev) +
  coord_flip()

plot_unimodal

# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/Figures/output/unimodal_performance.pdf",
#        plot_unimodal,
#        width = 23, height = 15, units = "cm")


df_combi <- data.frame(Data = factor(c("clinical + \nPRS", "clinical + \ntranscriptomic", "clinical + \nproteomic", "clinical + \nlipids and \nmetabolites", "all modalities"), 
                                     levels = c("clinical + \nPRS", "clinical + \ntranscriptomic", "clinical + \nproteomic", "clinical + \nlipids and \nmetabolites", "all modalities")), 
                       Accuracy = c(0.61, 0.64, 0.72, 0.66, 0.67),
                       AUROC = c(0.64, 0.61, 0.78, 0.65, 0.70))

plot_combi <- ggplot(df_combi, aes(Data, AUROC)) + 
  geom_point(aes(colour=AUROC), size = 5) +
  scale_color_continuous(low = "yellow", high = "red", na.value = NA, limits=c(0.5, 0.8)) +
  theme_light() +
  ylim(0.5, 1.0) +
  geom_text(aes(label = format(AUROC, nsmall = 2), vjust = AUROC-2), color = "black", size=6) +
  theme(text = element_text(size=24),
        axis.text=element_text(colour="black"), 
        legend.justification = 'top') +
  scale_x_discrete(limits=rev) +
  coord_flip()

plot_combi

# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/Figures/output/multimodal_performance.pdf",
#        plot_combi,
#        width = 23, height = 15, units = "cm")



  

