### Model Eval Results
mod_eval_results = mer
mer = mod_eval_results

### Subsetting model structures 6-11
nll_eval_results = mod_eval_results[, 2:7]

## Determining the number of NA iterations
na_index = c()
for(i in 1:length(nll_eval_results[ ,1])){
  if(any(is.na(nll_eval_results[i, ]))){
    na_index = c(na_index, i)
  }
}
na_index = unique(na_index)
### How many iterations failed to converge?
length(na_index)
 # nll_eval_results = nll_eval_results[-na_index, ]

## Getting summary stats for NLL models
range(nll_eval_results, na.rm = TRUE)
nll_vec = as.vector(nll_eval_results)
nll_vec = nll_vec[!is.na(nll_vec)]
mean(nll_vec)
sd(nll_vec)

### Determining which model performed best over all iterations
best_models = c()
for(i in 1:dim(nll_eval_results)[1]){
  best_models = c(best_models, names(which.min(nll_eval_results[i, ])))
}
table(best_models)

### Getting stats on the best performing model
best_nll_mod = mod_eval_results[ ,best_model]
range(best_nll_mod, na.rm = TRUE)
mean(best_nll_mod, na.rm = TRUE)
sd(best_nll_mod, na.rm = TRUE)

### Getting stats on the model based only on tagging data
mod_5 = as.vector(mer[ ,'model 5'])
range(mod_5, na.rm = TRUE)
mean(mod_5[!is.na(mod_5)])
sd(mod_5[!is.na(mod_5)])

### Comparing the perfered model to the tagging data only model
tagging_vs_composite_df = cbind(mer$`model 5`, mod_eval_results[ ,best_model])
colnames(tagging_vs_composite_df) = c('model 5', best_model)

tagging_vs_composite = c()
for(i in 1:length(tagging_vs_composite_df[ ,1])){
  tagging_vs_composite = c(tagging_vs_composite, colnames(tagging_vs_composite_df)[which.min(tagging_vs_composite_df[i, ])])
}
table(tagging_vs_composite)

### Summary stats on tagging and integrative models
pred_var_diff_tvc = tagging_vs_composite_df[ ,1] - tagging_vs_composite_df[ ,2]
range(pred_var_diff_tvc, na.rm = TRUE)
mean(pred_var_diff_tvc, na.rm = TRUE)
sd( as.vector(pred_var_diff_tvc)[!is.na(as.vector(pred_var_diff_tvc))])

#### Getting summary stats on all literature models 
lit = mod_eval_results[, 8:18]
range(lit, na.rm = TRUE)
mean(lit, na.rm = TRUE)
lit_vec = as.vector(lit)
lit_vec = lit_vec[!is.na(lit_vec)]
mean(lit_vec)
sd(lit_vec)

## Comparing Literatuere, MLE, and Bayesian models
model_structure_selection = data.frame()
nll_names = colnames(pref_mod)
lit_names = colnames(mod_eval_results)[8:18]
bayes_names = colnames(mod_eval_results)[19:22]
for(i in 1:length(mod_eval_results[ ,1])){
  score_ens = min(pref_mod[i], na.rm = TRUE)
  best_ens = best_model
  score_lit = min(mod_eval_results[i,8:18], na.rm = TRUE)
  best_lit = lit_names[which(mod_eval_results[i,8:18] == score_lit)]
  score_bayes = min(mod_eval_results[i,19:22], na.rm = TRUE)
  best_bayes = bayes_names[which(mod_eval_results[i,19:22] == score_bayes)]
  best_overall = c('MLE', 'Lit', 'Bayes')[which.min(c(score_ens, score_lit, score_bayes))]
  best_mod = c(best_ens, best_lit, best_bayes)[which.min(c(score_ens, score_lit, score_bayes))]
  write_line = data.frame('best_ll_mod' = best_ens, 'score_ensemble' = score_ens, 'best_lit_mod' = best_lit, 'score_lit' = score_lit, 'best_bayes_mod' = best_bayes, 'score_bayes' = score_bayes, 'best_model' = best_overall, 'best_mod' = best_mod)
  model_structure_selection = rbind(model_structure_selection, write_line)
}
lit_eval_results_table = aggregate(model_structure_selection$best_lit_mod, by = list(model_structure_selection$best_lit_mod), FUN = length)
best_lit_mod = lit_eval_results_table$Group.1[which.max(lit_eval_results_table$x)]

### Getting summary stats on the best performing literature model
best_lit = mod_eval_results[ ,as.character(best_lit_mod)]
range(best_lit, na.rm = TRUE)
best_lit_vec = as.vector(best_lit)
best_lit_vec = best_lit_vec[!is.na(best_lit_vec)]
mean(best_lit_vec)
sd(best_lit_vec)

## Getting summary stats for the second best performing literature model
second_best_lit_mod = as.character(lit_eval_results_table$Group.1[order(lit_eval_results_table$x, decreasing = TRUE)[2]])
second_best_lit = mod_eval_results[ ,as.character(second_best_lit_mod)]
second_best_lit_vec = as.vector(second_best_lit)
range(second_best_lit_vec, na.rm = TRUE)
mean(second_best_lit_vec, na.rm = TRUE)
sd(second_best_lit_vec, na.rm = TRUE)
