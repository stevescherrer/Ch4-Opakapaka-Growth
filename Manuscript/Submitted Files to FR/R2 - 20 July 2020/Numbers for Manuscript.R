### Getting required numbers for manuscript
load(file.path(results_dir, 'Paka VBGF Workspace For FR Revisions.RData'))

## 3.2 Estimating Growth Parameters from Tagging Data: Bayesian Approach

# Model 1, which incorporated individual variability in both L_∞ and K, yielded mean parameter estimates of...
round(model_1$BUGSoutput$mean$Linf_mu, digits = 2)
# 61.49
round(cv_df[cv_df$model == 'model_1' & cv_df$Parameter == "Linf", 'cv'], digits = 2)
#  2.56
round(model_1$BUGSoutput$mean$k_mu, digits = 2)
# 0.3
round(cv_df[cv_df$model == 'model_1' & cv_df$Parameter == "k", 'cv'], digits = 2)
# 8.33

# L_∞ and K parameter estimates for Model 2, where K was fixed, were ...
round(model_2$BUGSoutput$mean$Linf_mu, digits = 2)
# 61.61
round(cv_df[cv_df$model == 'model_2' & cv_df$Parameter == "Linf", 'cv'], digits = 2)
# 2.72
round(model_2$BUGSoutput$mean$k_mu, digits = 2)
# 0.29
round(cv_df[cv_df$model == 'model_2' & cv_df$Parameter == "k", 'cv'], digits = 2)
# 45.61

# Under Model 3, where L_∞ was fixed and K was fit freely ...
round(model_3$BUGSoutput$mean$Linf_mu, digits = 2)
# 74.14
round(cv_df[cv_df$model == 'model_3' & cv_df$Parameter == "Linf", 'cv'], digits = 2)
# 41.03
round(model_3$BUGSoutput$mean$k_mu, digits = 2)
# 0.17
round(cv_df[cv_df$model == 'model_3' & cv_df$Parameter == "k", 'cv'], digits = 2)
# 8.67

# ... for Model 4, where both parameters were fixed.
round(model_4$BUGSoutput$mean$Linf_mu, digits = 2)
# 73.62
round(cv_df[cv_df$model == 'model_4' & cv_df$Parameter == "Linf", 'cv'], digits = 2)
# 42.71
round(model_4$BUGSoutput$mean$k_mu, digits = 2)
# 0.18
round(cv_df[cv_df$model == 'model_4' & cv_df$Parameter == "k", 'cv'], digits = 2)
# 72.91

## DIC values
# Model 4 had the lowest DIC ...
model_4$BUGSoutput$DIC
  # 4788.539
# ...followed by Model 3...
model_3$BUGSoutput$DIC
  # 4974.611
# ... , and Model 2 ...
model_2$BUGSoutput$DIC
  # 8518.045
# ... while Model 1 had the highest DIC
model_1$BUGSoutput$DIC
  # 8581.735

### 3.4 Comparing model performance

# Across all 10,000 cross validation iterations to determine the preffered integrative model structure, the six candidate models produced RMSE values that ranged between... 
round(range(c(mod_eval_results$`model 6`, mod_eval_results$`model 7`, mod_eval_results$`model 8`, mod_eval_results$`model 9`, mod_eval_results$`model 10`, mod_eval_results$`model 11`), na.rm = T), 2)
  #  2.78 4.95
round(mean(c(mod_eval_results$`model 6`, mod_eval_results$`model 7`, mod_eval_results$`model 8`, mod_eval_results$`model 9`, mod_eval_results$`model 10`, mod_eval_results$`model 11`), na.rm = T), 2)
  # 3.9
round(sd(c(mod_eval_results$`model 6`, mod_eval_results$`model 7`, mod_eval_results$`model 8`, mod_eval_results$`model 9`, mod_eval_results$`model 10`, mod_eval_results$`model 11`), na.rm = T), 2)
  # 0.3

# The structure of Model 11 outperformed competing models during cross validation...
length(which(integrative_model_scores == 'model 11'))
  # 2192

# RMSE for this model ranged between...
round(range(mod_eval_results$`model 11`, na.rm = T), digits = 2)
 # 2.86 4.89
round(mean(mod_eval_results$`model 11`, na.rm = T), digits = 2)
  # 3.89
round(sd(mod_eval_results$`model 11`, na.rm = T), digits = 2)
  # 0.29

# The structure of Model 11 performed better than the structure of Model 5 during cross validation... 
length(which(overall_max_likelihood_models == 'model 11'))
  # 5672

# Differences in RMSE between these two competing structures ranged between 
round(range((mod_eval_results$`model 11` - mod_eval_results$`model 5`), na.rm = T), 2) 
  # -1.15  0.14
round(mean(range((mod_eval_results$`model 11` - mod_eval_results$`model 5`), na.rm = T)), 2)
  # -0.51
round(sd((mod_eval_results$`model 11` - mod_eval_results$`model 5`), na.rm = T),2)
  # 0.11

# with structure of Model 5, fit exclusively using tagging data, producing RMSE values that ranged between...
round(range(mod_eval_results$`model 5`, na.rm = T), 2)
  # 2.84 5.28
round(mean(mod_eval_results$`model 5`, na.rm = T), 2)
  # 3.9
round(sd(mod_eval_results$`model 5`, na.rm = T), 2)
  # 0.33

# 3.5 Sensitivity Analysis
# Parameters estimated using the observed and synthetic data differed by as much as ...
round(max(sensitivity_table$percent_k, sensitivity_table$percent_linf),2)
  # 205.08
round(median(c(abs(sensitivity_table$percent_k), abs(sensitivity_table$percent_linf))),2)
  # 21.28

## Model 7 was the most robust, with the smallest difference between L_∞ and K parameters estimated with real and synthetic datasets 
sensitivity_table[sensitivity_table$Model == 'Model 7', 'percent_linf']
# -1.13
sensitivity_table[sensitivity_table$Model == 'Model 7', 'percent_k']
# 8.7

## The preferred integrated model (Model 11) was the second best performing model overall with L_∞ differing by 2.45% and K by 15.47%. 
sensitivity_table[sensitivity_table$Model == 'Model 11', 'percent_linf']
  # -2.45
sensitivity_table[sensitivity_table$Model == 'Model 11', 'percent_k']
  # 15.47

# Parameter estimates for Model 1, the Bayesian model that accounted for individual differences in each parameter and had the lowest coefficient of variation across both parameters, differed by...
abs(sensitivity_table[sensitivity_table$Model == 'Bayesian Model 1 (with Priors)', 'percent_linf'])
  # 36.8
abs(sensitivity_table[sensitivity_table$Model == 'Bayesian Model 1 (with Priors)', 'percent_k'])
  # 122.46

## Parameters for Model 4, the Bayesian model with the lowest DIC score, differed between observed and synthetic data by ...
abs(sensitivity_table[sensitivity_table$Model == 'Bayesian Model 4 (with Priors)', 'percent_linf'])
  # 20.6
abs(sensitivity_table[sensitivity_table$Model == 'Bayesian Model 4 (with Priors)', 'percent_k'])
  # 135.04

### Table 3 values to plug and play
table_3 = data.frame(
  rbind(
    cbind('Model 1', model_1$BUGSoutput$summary[ ,colnames(model_1$BUGSoutput$summary) %in% c('mean', 'sd', '2.5%', '50%', '97.5%', 'Rhat', 'n.eff')]),
    cbind('Model 2', model_2$BUGSoutput$summary[ ,colnames(model_1$BUGSoutput$summary) %in% c('mean', 'sd', '2.5%', '50%', '97.5%', 'Rhat', 'n.eff')]),
    cbind('Model 3', model_3$BUGSoutput$summary[ ,colnames(model_1$BUGSoutput$summary) %in% c('mean', 'sd', '2.5%', '50%', '97.5%', 'Rhat', 'n.eff')]),
    cbind('Model 4', model_4$BUGSoutput$summary[ ,colnames(model_1$BUGSoutput$summary) %in% c('mean', 'sd', '2.5%', '50%', '97.5%', 'Rhat', 'n.eff')])
  )
)
View(table_3)

## Write out for easy drag & drop
write.csv(table_3, file.path(results_dir, 'Tables/table_3 - Bayesian Parameter Estimates.csv'))

