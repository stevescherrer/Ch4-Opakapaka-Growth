###### Generating Images For Growth Paper
proj_dir = '/Volumes/GoogleDrive/My Drive/Weng Lab/Personal_Folders/Steve/dissertation work/Ch 4. Opakapaka Growth/Analysis'
results_dir = file.path(proj_dir, 'results')
lit_vbgc_params = read.csv(file.path(proj_dir, 'data/lit_vbgf_params.csv'), stringsAsFactors = FALSE)
lit_vbgc_params = lit_vbgc_params[!is.na(lit_vbgc_params$Linf), ]
colnames(lit_vbgc_params) = c('author', 'n', 'linf', 'k', 't0', 'region', 'method')


#### Some general Functions
std_error = function(x){
  #### Calculates standard error of set (x)
  sqrt(var(x)/length(x))
}

von_b_eq = function(t, t0, linf, k){
  ## Get estimated length at time t using von bertalanffy function
  return(linf*(1-exp(-k*(t-t0))))
}

predict_recapture_length = function(Lm, dt, linf = 65.95546, k = 0.2369113){
  ## Get estimated length at recapture of a given individual using von Bertalanffy function as paramterized by Faben
  return(Lm + (linf - Lm) * (1 - exp(-k * dt)))
}

plot_error_bars = function(x, y = NULL, mean_y = NULL, se_y = NULL, color = 'black'){
  if(is.null(y) & is.null(mean_y) & is.null(se_y)){
    print('Either argument "y" or arguments "mean_y" and "se_y" must not be NULL')
  }else if (is.null(y)){
    arrows(x, mean_y - (se_y), x, mean_y + (se_y), length=0.05, angle=90, code=3, col = color, lwd = 1)
  }else if (!is.null(y)){
    arrows(x, mean(y) - se(y), x, mean(y) + se(y), length=0.05, angle=90, code=3, col = color, lwd = 1)
  }
}


##### Plotting Figures ####

#### Histogram of length at tagging and recapture 
pdf(file.path(results_dir, 'Figure 1: Tag recapture length histogram.pdf'))
par(mfrow = c(2, 1))
hist(paka_growth$Lm, col = rgb(1,0,0,0.5), xlim = c(0, 100), ylim = c(0, 200), xlab = 'Reported Fork Length (cm)', ylab = "n Individuals", main = NULL)
hist(paka_growth$Lr, col = rgb(0,0,1,0.5), add = TRUE)
legend('topright', c("Tagged", "Recaptured"), fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))
hist(paka_growth$dt, col = 'grey', xlim = c(0, 12), ylim = c(0, 100), xlab = 'Time at Liberty (years)', ylab = "n Individuals", main = NULL, breaks = seq(0, 11, .5))
dev.off()


#### Plotting predicted vs. observed length at recapture 
#### For our Modesl

### Constructing data frame of predicted values and residuals under each model 
predicted_recapture_lengths = data.frame(NULL)
predicted_recapture_residuals = data.frame(NULL)

models_to_predict = lit_vbgc_params[lit_vbgc_params$author %in% c(paste("Bayesian Model", 1:2),  "Maximum Likelihood - Model 5", "Maximum Likelihood - Integrative Model"), ]
models_to_predict$linf = as.numeric(models_to_predict$linf)
models_to_predict$k = as.numeric(models_to_predict$k)
models_to_predict$t0 = as.numeric(models_to_predict$t0)
  models_to_predict$t0[is.na(models_to_predict$t0)] = 0

for(i in 1:nrow(models_to_predict)){
  predicted_recapture_lengths = rbind(predicted_recapture_lengths, predict_recapture_length(Lm = tagdat[ ,1], dt = tagdat[ ,4], linf = as.numeric(models_to_predict$linf[i]), k = as.numeric(models_to_predict$k[i])))
  predicted_recapture_residuals = rbind(predicted_recapture_residuals, predict_recapture_length(Lm = tagdat[ ,1], dt = tagdat[ ,4], linf = as.numeric(models_to_predict$linf[i]), k = as.numeric(models_to_predict$k[i])) - tagdat[ ,2])
}

rownames(predicted_recapture_lengths)   = c('Model 1', 'Model 2', 'Model 5', 'Model 11')
rownames(predicted_recapture_residuals) = c('Model 1', 'Model 2', 'Model 5', 'Model 11')

pdf(file.path(results_dir, 'Figure 3 - Predicted vs. Observed LR with validation data.pdf'), height = 8.5, width = 8.5)
par(mfrow = c(2, 2))
for(i in 1:nrow(predicted_recapture_lengths)){
  model_id = row.names(predicted_recapture_lengths)[i]
  plot(y = predicted_recapture_lengths[i, ], x = tagdat[ ,2],
       xlab = 'Observed Recapture FL (cm)', xlim = c(15, 80), 
       ylab = 'Predicted Recapture FL (cm)', ylim = c(15, 80),
       main = paste(model_id,'\n Linf = ', models_to_predict$linf[i],', K = ', models_to_predict$k[i], sep = ""))
  abline(1, 1, lty = 2)
  arrows(x0 = 60, y0 = 70, x1 = 65, y1 = 70, length = 0.1, angle = 30, code = 2)
  text(x = 45, y = 70, labels = "Line of 1:1 agreement", cex = .75)
  model_var = sum((predicted_recapture_lengths[i, ] - tagdat_validate[ ,2])^2) / length(tagdat_validate[ ,2])
  # text(x = 60, y = 30, labels = paste("Predictive Variance:", round(model_var, digits = 3)), cex = .75)
}
dev.off()


#### Plot of VBGF Curves
### For this study
pdf(file.path(results_dir, 'Figure 4 - VBGF Plots for Bayesian and MLE models.pdf'), height = 8.5, width = 11)
par(mfrow = c(1, 1))

## Define how far the X axis extends. 
max_age = 20

## Set up an empty plot
plot(x = c(0, max_age), y = c(0, 100),
     xlim = c(0, max_age), xlab = 'Years',
     ylim = c(0, 75), ylab = "Fork Length (cm)",
     cex = 1.5, cex.axis = 1.5, cex.lab = 1.5, 
     type = 'n')

## Add extra tick marks
axis(1, at = seq(0, 20, 2.5), labels = F)
axis(2, at = seq(0, 70, 5), labels = F)

## Loop through each model, draw a predicted line for that model on the plot
plot_col = rainbow(length(models_to_predict$author))
for(i in 1:length(models_to_predict$author)){
  estimates = von_b_eq(t = seq(0, max_age, .25), t0 = 0, linf = models_to_predict$linf[i], k = models_to_predict$k[i])
  lines(estimates ~ seq(0, max_age, .25), col = plot_col[i], lwd = 3)
}

## Adding lines representing smallest and largest fish tagged at capture
abline(min(tagdat[ ,1]), 0, lty = 2)
abline(max(tagdat[ ,1]), 0, lty = 2)
text(x = 0.5, y = min(tagdat[ ,1]) + 5, labels = 'min Lm', cex = 1.5)
text(x = 0.5, y = max(tagdat[ ,1]) + 5, labels = 'max Lm', cex = 1.5)

## Adding a legend showing which line corrosponds to which model
legend('bottomright', legend = row.names(predicted_recapture_lengths), lty = c(1, 1, 1, 1, 1, 1, 1, 1, 2), col = c(plot_col, 'black', 'black'), cex = 1.5, lwd = 3)
dev.off()

### For all studies
# pdf(file.path(results_dir, 'Fig. 4 - VBGF Plots for all studies.pdf'), height = 8.5, width = 11)
# plot(x = c(0, 46), y = c(0, 100), type = "n",
#      xlim = c(0, 46), xlab = 'Age (years)',
#      ylim = c(0, 100), ylab = "Length (cm)",
#      main = "VBGF for P. Filamentosus", cex = 1.5, cex.axis = 1.5, cex.lab = 1.5)
# 
# plot_col = rainbow(length(lit_vbgc_params$author))
# for(i in 1:length(lit_vbgc_params$author)){
#   estimates = von_b_eq(0:46, t0 = lit_vbgc_params$t0[i], linf = lit_vbgc_params$linf[i], k = lit_vbgc_params$k[i])
#   lines(estimates ~ c(0:46), col = plot_col[i])
# }
# 
# legend('bottomright', legend = lit_vbgc_params$author, lty = c(1, 1, 1, 1, 1, 1, 1, 1, 2), col = c(plot_col, 'black', 'black'), cex = 1)
# dev.off()


#### Plotting growth vs time at liberty
pdf(file.path(results_dir, 'delta_t vs delta_L plot.pdf'), height = 8.5, width = 11)
plot(x = paka_growth$dt, y = paka_growth$dL,
     ylab = 'ΔL (cm)', ylim = c(0, 50),
     xlab = 'ΔT (years)', xlim = c(0, 11),
     main = 'Comparing Individual Growth and Time at Liberty',
     pch = 19)
abline(lm(paka_growth$dL ~ paka_growth$dt), col = 'red')
lm_sum = summary(lm(paka_growth$dL ~ paka_growth$dt))
text(x = 8, y = 10, labels = paste('R^2 =', round(lm_sum$r.squared, digits = 2)))
dev.off()

#### Histogram of time at liberty
pdf(file.path(results_dir, 'delta_t and delta_L histograms.pdf'), height = 8.5, width = 11)
par = mfrow(c(2, 1))
hist(as.numeric(paka_growth$dt), ylim = c(0, 100), xlab = 'Δt (years)', main = NULL, right = TRUE, col = 'Grey', breaks = seq(from = 0, to = 12, by = .5))
#### Histogram of growth at liberty
hist(paka_growth$dL, main = NULL, xlab = 'ΔL (cm)', xlim = c(-1, 45), ylim = c(0, 35), right = TRUE, col = 'Grey', breaks = seq(from = -1, to = 50, by = 1))
dev.off()




#################### OLD PLOTS ###########################


# #### Plotting residual vs. observed residuals for each model
# pdf(file.path(results_dir, 'Figure 4 - Residual vs. Observed LR with validation data.pdf'), height = 8.5, width = 11)
# par(mfrow = c(2, 2))
# for(i in 1:nrow(predicted_recapture_lengths)){
#   model_id = row.names(predicted_recapture_lengths)[i]
#   plot(as.numeric(predicted_recapture_residuals[i, ]) ~ tagdat[ ,2],
#        xlab = 'Observed FL (cm)', xlim = c(15, 80),
#        ylab = 'Model Residual (cm)', ylim = c(-1*max(abs(range(predicted_recapture_residuals))), max(abs(range(predicted_recapture_residuals)))),
#        main = paste(model_id,'\n Linf = ', models_to_predict$linf[i],', k = ', models_to_predict$k[i], sep = ""))
#   abline(0, 0, lty = 2)
# }
# dev.off()

#### Plotting predicted vs. observed for all studies
lsandks = data.frame()

predicted_recapture_lengths = data.frame(NULL)
predicted_recapture_residuals = data.frame()
for(i in 1:length(lit_vbgc_params$author)){
  med_boot_l = lit_vbgc_params$linf[i]
  med_boot_k = lit_vbgc_params$k[i]
  med_boot_a0 = lit_vbgc_params$t0[i]
  predicted_recapture_lengths = rbind(predicted_recapture_lengths, predict_recapture_length(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = med_boot_l, k = med_boot_k))
  predicted_recapture_residuals = rbind(predicted_recapture_residuals, predict_recapture_length(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = med_boot_l, k = med_boot_k) - tagdat_validate[ ,2])
}
rownames(predicted_recapture_lengths)   = lit_vbgc_params$author
rownames(predicted_recapture_residuals) = lit_vbgc_params$author

pdf(file.path(results_dir, 'Predicted vs. Observed LR - all studdies.pdf'), height = 8.5, width = 11)
par(mfrow = c(ceiling(length(rownames(lit_vbgc_params))/4), 6))
for(i in 1:nrow(predicted_recapture_lengths)){
  plot(y = predicted_recapture_lengths[i, ], x= tagdat_validate[ ,2],
       xlab = 'Observed Length at Recapture (cm)', xlim = c(15, 80), 
       ylab = 'Predicted Length at Recapture (cm)', ylim = c(15, 80),
       main = paste('Linf = ', round(lit_vbgc_params$linf[i], digits = 2),'k = ', round(lit_vbgc_params$k[i], digits = 2), sep = ""))
  abline(1, 1, lty = 2)
  arrows(x0 = 60, y0 = 70, x1 = 65, y1 = 70, length = 0.1, angle = 30, code = 2)
  text(x = 45, y = 70, labels = "Line of 1:1 agreement", cex = .75)
  text(x = 50, y = 30, labels = paste("Correlation Coefficient:", round(var(x = as.numeric(predicted_recapture_lengths[i, ]), y = as.numeric(tagdat_validate[ ,2])), digits = 3)), cex = .75)
}
dev.off()

### Plotting predicted vs. observed residuals - All Studies
pdf(file.path(results_dir, 'Residual vs. Observed LR - all studies.pdf'), height = 8.5, width = 11)
par(mfrow = c(ceiling(length(rownames(lit_vbgc_params))/4), 6))
for(i in 1:nrow(predicted_recapture_residuals)){
  plot(as.numeric(predicted_recapture_residuals[i, ]) ~ tagdat_validate[ ,2],
       xlab = 'Observed Length at Recapture (cm)', xlim = c(15, 80),
       ylab = 'Residual Length at Recapture (cm)', ylim = c(-1*max(abs(range(predicted_recapture_residuals))), max(abs(range(predicted_recapture_residuals)))),
       main = paste('Linf = ', round(lit_vbgc_params$linf[i], digits = 2),'k = ', round(lit_vbgc_params$k[i], digits = 2), sep = ""))
  abline(0, 0, lty = 2)
}
dev.off()



#### CAN WE DO ANOVA?
pds = data.frame(NULL)
for(i in 1:length(names(bootstrap_results))){
  new_df = data.frame(bootstrap_results[[i]]$raw_boot_data$mu.L)
  new_df$model = as.factor(names(bootstrap_results)[i])
  pds = rbind(pds, new_df)
}

param_anov = aov(pds$bootstrap_results..i...raw_boot_data.mu.L~pds$model)
TukeyHSD(param_anov)
plot(param_anov)





#### Analysis using Patterson III et al 2001 as a guide 


#### Making a plot of predicted total length compared to observed total length at recatpure
#### Constructing a vector of daily predicted growth using vonB params, then comparing tagging results
pdf(file.path(results_dir, 'predicted vs. observed recapture lengths.pdf'), height = 8.5, width = 11)
par(mfrow = c(2, 3))

## Subsetting data to only our models
subset_data = lit_vbgc_params[lit_vbgc_params$author %in% c('Bayesian Model 1', 'Bayesian Model 2', 'Bayesian Model 3', 'Bayesian Model 4', 'Maximum Likelihood - Tagging Data', 'Maximum Likelihood - Tagging & Otolith Data'), ]
for(i in 1:length(subset_data$author)){
  print(i)
  ## Selecting indvidual model parameter estimates
  model_id = subset_data$author[i]
  linf = subset_data$linf[i]
  k = subset_data$k[i]

  ## Predicting length at recapture
  paka_growth$Lr_predicted = predict_recapture_length(Lm = paka_growth$Lm, dt = paka_growth$dt, linf = linf, k = k)

  ## Modeling predicted and observed growth
  pred_obs_lm = lm(paka_growth$Lr_predicted~paka_growth$Lr)
  # summary(pred_obs_lm)
  # plot(pred_obs_lm)
  plot(paka_growth$Lr_predicted ~ paka_growth$Lr, 
       main = model_id,
       ylab = 'Predicted FL at Recapture (cm)', xlab = 'Observed FL at Recapture (cm)',
       xlim = c(15, 80), ylim = c(15, 80), pch = 19)
  abline(1, 1, lty = 2)
  arrows(x0 = 60, y0 = 70, x1 = 65, y1 = 70, length = 0.1, angle = 30, code = 2)
  text(x = 43, y = 70, labels = "Line of 1:1 agreement")
  text(x = 60, y = 20, labels = paste('R^2: ', round(summary(pred_obs_lm)$r.squared, digits = 3)))
}
dev.off()





## Modeling and Plotting our predictions vs. those from previous studies in the MHI

predicted_recapture_lengths = list()
for(i in 1:length(results$Data)){
  predicted_recapture_lengths[[i]] = predict_recapture_length(Lm = paka_growth$Lm, dt = paka_growth$dt, linf = lit_vbgc_params$linf[i], k = lit_vbgc_params$k[i])
  names(predicted_recapture_lengths)[[i]] = as.character(lit_vbgc_params$author[i])
}

pdf(file.path(results_dir, 'Comparision to Literature Estimates.pdf'), height = 11, width = 8.5)
par(mfrow = c((length(lit_vbgc_params$author) / 2), 2), oma = c(0,0,2,0))
for(i in 1:(length(predicted_recapture_lengths)-1)){
  ## Making Plot
  plot(predicted_recapture_lengths$`Scherrer, Franklin 2017` ~ predicted_recapture_lengths[[i]], 
       ylab = 'Scherrer and Franklin', xlab = lit_vbgc_params$author[i], pch = 19,
       xlim = c(15, 75), ylim = c(15, 75))  
  abline(1, 1, lty = 2)
  ## Modeling
  pred_lm = lm(predicted_recapture_lengths$`Scherrer, Franklin 2017` ~ predicted_recapture_lengths[[i]])
  text(x = 60, y = 25, labels = c(paste('Slope: ', round(summary(pred_lm)$coefficients[2, 1], digits = 3), "+/-", round(summary(pred_lm)$coefficients[2, 2], digits = 3))))
  text(x = 60, y = 20, labels = c(paste('R^2: ', round(summary(pred_lm)$r.squared, digits = 3))))
  print(cor(predicted_recapture_lengths$`Scherrer, Franklin 2017`, predicted_recapture_lengths[[i]]))
}
mtext('Comparing Predicted Fork Length at Recapture (cm)', outer = TRUE, side = 3)
dev.off()









################# SOME OTHER SHIT ##################






#### Making a table to compare Rsqs.
predicted_recapture_lengths = data.frame(NULL, stringsAsFactors = F)
predicted_recapture_residuals = data.frame(NULL, stringsAsFactors = F)

for(i in 1:length(names(bootstrap_results))){
  med_boot_l = bootstrap_results[[i]]$boot_stats['Median','mu.L']
  med_boot_k = bootstrap_results[[i]]$boot_stats['Median','k']
  med_boot_a0 = bootstrap_results[[i]]$boot_stats['Median','mu.A']
  predicted_recapture_lengths = rbind(predicted_recapture_lengths, predict_recapture_length(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = bootstrap_results[[i]]$boot_stats['Median','mu.L'], k = bootstrap_results[[i]]$boot_stats['Median','k']))
  predicted_recapture_residuals = rbind(predicted_recapture_residuals, predict_recapture_length(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = bootstrap_results[[i]]$boot_stats['Median','mu.L'], k = bootstrap_results[[i]]$boot_stats['Median','k']) - tagdat_validate[ ,2])
}

### Adding in literature models
rsq_compare_df = data.frame(stringsAsFactors = F)
var_compare_df = data.frame(stringsAsFactors = F)
for(i in 1:length(lit_vbgc_params$author)){
  predicted_recapture_lengths = rbind(predicted_recapture_lengths, predict_recapture_length(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = lit_vbgc_params$linf[i], k = lit_vbgc_params$k[i]))
  predicted_recapture_residuals = rbind(predicted_recapture_residuals, predict_recapture_length(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = lit_vbgc_params$linf[i], k = lit_vbgc_params$k[i]))
}

rownames(predicted_recapture_lengths)   = c(lit_vbgc_params$author)
rownames(predicted_recapture_residuals) = c( lit_vbgc_params$author)

# rownames(predicted_recapture_lengths)   = c(names(bootstrap_results), lit_vbgc_params$author)
# rownames(predicted_recapture_residuals) = c(names(bootstrap_results), lit_vbgc_params$author)

for(i in 1:length(rownames(predicted_recapture_lengths))){
  model_id = as.character(rownames(predicted_recapture_lengths)[i])
  model_rsq = as.numeric(summary(lm(as.numeric(predicted_recapture_lengths[i, ]) ~ tagdat_validate[,2]))$adj.r.squared)
  model_var =  as.numeric(sum((as.numeric(predicted_recapture_lengths[i, ]) - tagdat_validate[,2])^2) / length(tagdat_validate[,2]))
  newline = data.frame('model id' = model_id, 'adj_r.sq' = model_rsq, 'variance' = model_var)
  print(newline)
  var_compare_df = rbind(var_compare_df, newline)
}
write.csv(var_compare_df, file.path(results_dir, 'model performance comparison table.csv'))

mean(var_compare_df$variance[2:7])
std_error(var_compare_df$variance[2:7])

range(var_compare_df$variance[7:26])
mean(var_compare_df$variance[7:26])
std_error(var_compare_df$variance[7:26])

############ FIGURES FOR TALKS ####################
load('/Users/stephenscherrer/Google Drive/Weng Lab/Data/Bottomfish/Okamoto_Mark_Recapture/results/FINAL 1000 run 2018-04-15 21:40:13/workspace_image_preboot.RData')

### Comparing 10000 iterations between tagging only and best integrative model
integrative_vs_tagging_model_scores_table = table(integrative_vs_tagging_model_scores)
integrative_vs_tagging_model_scores_table = integrative_vs_tagging_model_scores_table[2:1]
 pdf(file.path(results_dir, 'Barplot of Integrative vs. Tagging Models.pdf'), height = 8.5, width = 11)
  barplot(integrative_vs_tagging_model_scores_table, ylim = c(0, 7500), col = 'lightblue')
 dev.off()
 
 
### Comparing 10000 iterations between integrative models
 int_mod_results = mod_eval_results[, 2:7]
 best_int_mod = c()
 for(i in 1:length(int_mod_results[,1])){
   best_int_mod = c(best_int_mod, names(which.min(int_mod_results[i, ])))
 }
 int_mod_table = table(best_int_mod)
 table_order = c(3, 4, 5, 6, 1, 2)
 int_mod_table = int_mod_table[table_order]
 
 pdf(file.path(results_dir, 'Barplot of Integrative Model Performance.pdf'), height = 8.5, width = 11)
  barplot(height = int_mod_table, col = 'lightblue', ylim = c(0, 5000))
 dev.off()
