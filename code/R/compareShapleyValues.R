library(tidyverse)
library(reshape2)
setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/")

xgb_R_Shap1 <- read_csv("Predictions_explorative/R.scripts/data/analyte_importance_XGB-Shap.csv")

xgb_R_Shap2 <- read_csv("Predictions_explorative/R.scripts/data/analyte_shap_xgb.csv")

xgb_Python_Shap <- read_csv("Predictions_explorative/Python.scripts/output/avg_shap_abs.csv")
xgb_Python_Shap <- xgb_Python_Shap[,-1]

rank_table <- xgb_R_Shap1 %>% 
  left_join(xgb_R_Shap2, by = c("symbol","Analyte full name")) %>% 
  left_join(xgb_Python_Shap, by = c("symbol"="Protein"))

names(rank_table) <- c("protein", "protein_full_name", "shap_xgboost_v1", "shap_xgboost_v2", "shap_xgboost_python")

rank_table[is.na(rank_table$shap_xgboost_python),] #proteins left out by preprocessing have '0' value in other shap-explainers, so replace with 0
rank_table$shap_xgboost_python <- replace_na(rank_table$shap_xgboost_python, 0)

# Write to file
# write_csv(rank_table, "Predictions_explorative/R.scripts/2-yearChronicity/output/ShapRankTableXGBmodels.csv")

#plot relationships
plot(rank_table$shap_xgboost_v1, rank_table$shap_xgboost_v2, main = "Early model vs last best model in R",
     xlab = "shapleyValues_xgboost_v1", ylab = "shapleyValues_xgboost_v2",
     pch = 19, frame = FALSE)
abline(lm(shap_xgboost_v2 ~ shap_xgboost_v1, data = rank_table), col = "blue")

plot(rank_table$shap_xgboost_v2, rank_table$shap_xgboost_python, main = "R vs Python",
     xlab = "shapleyValues_xgboost_v2", ylab = "shapleyValues_xgboost_python",
     pch = 19, frame = FALSE)
abline(lm(shap_xgboost_python ~ shap_xgboost_v2, data = rank_table), col = "blue")

#correlations
cor.test(rank_table$shap_xgboost_v1, rank_table$shap_xgboost_v2, alternative = "two.sided", exact = FALSE, method = "pearson")
cor.test(rank_table$shap_xgboost_v2, rank_table$shap_xgboost_python, alternative = "two.sided", exact = FALSE, method = "pearson")

cormat <- cor(rank_table[,-c(1:2)], method = "pearson")
melted_cormat <- melt(cormat)
names(melted_cormat) <- c("explainers_x", "explainers_y", "correlation")
ggplot(data = melted_cormat, aes(x=explainers_x, y=explainers_y, fill=correlation)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  geom_text(aes(explainers_y, explainers_x, label = round(correlation,2)), color = "black", size = 4)











  