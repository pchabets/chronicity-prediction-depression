library(tidyverse)
library(pROC)
setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/")

# Load in models (ANN was done in python, values filled in manually in dataframe below)
eln_clin <- readRDS("eln_roc_2ychr_clin.RDS")
eln_prot <- readRDS("roc_2ychronicity_proteomics_elnet.RDS")
eln_prot_clin <- readRDS("roc_2ychronicity_prot_clin_ELNET.RDS")

svm_clin <- readRDS("svmL_roc_2ychr_clin.RDS")
svm_prot <- readRDS("roc_2ychronicity_proteomics_SVM-L.RDS")
svm_prot_clin <- readRDS("roc_2ychronicity_prot_clin_SVM-L.RDS")

rf_clin <- readRDS("rf_roc_2ychr_clin.RDS")
rf_prot <- readRDS("roc_2ychronicity_proteomics_RF.RDS")
rf_prot_clin <- readRDS("roc_2ychronicity_prot_clin_RF.RDS")

xgb_clin <- readRDS("roc_2ychronicity_clinical.RDS")
xgb_prot <- readRDS("roc_2ychronicity_proteomics.RDS")
xgb_prot_clin <- readRDS("roc_2ychronicity_prot_clin.RDS")

# Combine data in df
compare_df <- data.frame(model = factor(c(rep("ElasticNet", 3), 
                                   rep("SVM", 3),
                                   rep("RandomForest", 3),
                                   rep("XGBoost", 3),
                                   rep("ANN", 3)), levels = c("ElasticNet", "SVM", "RandomForest", "XGBoost", "ANN")),
                         
                         data = c(rep(c("clinical", "proteomics", "proteomics + clinical"), 5)),
                         
                         AUROC = c(auc(eln_clin), auc(eln_prot), auc(eln_prot_clin),
                                   auc(svm_clin), auc(svm_prot), auc(svm_prot_clin),
                                   auc(rf_clin), auc(rf_prot), auc(rf_prot_clin),
                                   auc(xgb_clin), auc(xgb_prot), auc(xgb_prot_clin),
                                   0.625, 0.652, 0.7315),
                         
                         CI_low = c(ci.auc(eln_clin, method = "delong")[1], ci.auc(eln_prot, method = "delong")[1], ci.auc(eln_prot_clin, method = "delong")[1],
                                    ci.auc(svm_clin, method = "delong")[1], ci.auc(svm_prot, method = "delong")[1], ci.auc(svm_prot_clin, method = "delong")[1],
                                    ci.auc(rf_clin, method = "delong")[1], ci.auc(rf_prot, method = "delong")[1], ci.auc(rf_prot_clin, method = "delong")[1],
                                    ci.auc(xgb_clin, method = "delong")[1], ci.auc(xgb_prot, method = "delong")[1], ci.auc(xgb_prot_clin, method = "delong")[1],
                                    0.53801123, 0.5529629, 0.64075974),
                         
                         CI_high = c(ci.auc(eln_clin, method = "delong")[3], ci.auc(eln_prot, method = "delong")[3], ci.auc(eln_prot_clin, method = "delong")[3],
                                     ci.auc(svm_clin, method = "delong")[3], ci.auc(svm_prot, method = "delong")[3], ci.auc(svm_prot_clin, method = "delong")[3],
                                     ci.auc(rf_clin, method = "delong")[3], ci.auc(rf_prot, method = "delong")[3], ci.auc(rf_prot_clin, method = "delong")[3],
                                     ci.auc(xgb_clin, method = "delong")[3], ci.auc(xgb_prot, method = "delong")[3], ci.auc(xgb_prot_clin, method = "delong")[3],
                                     0.71218411, 0.75048161, 0.82233485)
                         )


# Plot AUROC and CI for all models and data 
models_plot <- ggplot(compare_df, aes(x = data, y = AUROC, fill = model, colour = model)) +
  geom_point(position=position_dodge(width=0.5)) +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=.1, position=position_dodge(width=0.5)) +
  theme_gray(base_size = 24)
models_plot
# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/Figures/output/linear_non_linear_comparison_plot.pdf",
#        models_plot,
#        width = 35, height = 15, units = "cm")

models_plot2 <- ggplot(compare_df, aes(x=data, y=AUROC, group=model, color=model)) + 
  geom_line(position = position_dodge(width=0.5),
            linetype = 'dotted', 
            size = .75) +
  geom_point(position = position_dodge(width=0.5),
             size = 4)+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), 
                width=.25,
                size = .75,
                position=position_dodge(0.5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_gray(base_size = 24) +
  theme(legend.justification = 'top')
models_plot2
# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/Figures/output/linear_non_linear_comparison_plot2.pdf",
#        models_plot2,
#        width = 35, height = 15, units = "cm")

# Plot AUROC and CI for XGBoost models and data 
sub_df <- subset(compare_df, model == "XGBoost")
xgb_plot <- ggplot(sub_df, aes(x = data, y = AUROC, group = model, colour=model)) +
  geom_point(size=5)+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=.075, size=1) +
  geom_line(linetype = 'dotted', size=1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_gray(base_size = 24) + 
  expand_limits(y =c(0.5, 0.9)) +
  theme(legend.justification = 'top')
xgb_plot

# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/Figures/output/xgb_clin_prot_comparison.pdf",
#        xgb_plot,
#        width = 35, height = 15, units = "cm")

# ggplot(sub_df, aes(x=data, y = AUROC, colour = model)) +
#   geom_point(position=position_dodge(width=0.5)) +
#   # scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
#   geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=.1, position=position_dodge(width=0.5)) +
#   theme_gray(base_size = 16)




# Test significance of differences in AUROC
roc.test(eln_clin, eln_prot_clin, method = "delong", paired = FALSE)
roc.test(svm_clin, svm_prot_clin, method = "delong", paired = FALSE)
roc.test(rf_clin, rf_prot_clin, method = "delong", paired = FALSE)
roc.test(xgb_clin, xgb_prot_clin, method = "delong", paired = FALSE)

roc.test(xgb_prot_clin, eln_prot_clin, method = "delong", paired = TRUE)

# Test for extended clinical data
xgb_clin_ext <- readRDS("xgb_clinical_extended.RDS")
xgb_clin_ext_dinga <- readRDS("roc_2ychr_dinga_extended_clinical.RDS")
roc.test(xgb_clin, xgb_clin_ext, method = "delong", paired = FALSE)
roc.test(xgb_clin, xgb_clin_ext_dinga, method = "delong", paired = FALSE)
roc.test(xgb_clin_ext, xgb_clin_ext_dinga, method = "delong", paired = FALSE)
roc.test(xgb_prot_clin, xgb_clin_ext, method = "delong", paired = FALSE)
roc.test(xgb_prot_clin, xgb_clin_ext_dinga, method = "delong", paired = FALSE)

# Power calculations
power.roc.test(eln_prot_clin, rf_prot_clin)
power.roc.test(eln_prot_clin, xgb_prot_clin)

calc_power <- function(roc_clinical, roc_clin_prot, alpha = 0.05){
  # function to calculate power of alternative ROC under normal distribution
  critical_value <- qnorm(1-alpha/2, mean=auc(roc_clinical), sd=sqrt(var(roc_clinical)))
  power = 1 - pnorm(critical_value, mean=auc(roc_clin_prot), sd=sqrt(var(roc_clin_prot)))
  return(power)
}

calc_power(xgb_clin, xgb_prot_clin)
calc_power(rf_clin, rf_prot_clin)

crit_value_ann <- qnorm(1-0.05/2, mean=0.6250976715111737, sd=sqrt(0.0019742624291418723))
power_ann <- 1 - pnorm(crit_value_ann, mean=0.7315472936030617 , sd=sqrt(0.0021456380062029007))
power_ann




