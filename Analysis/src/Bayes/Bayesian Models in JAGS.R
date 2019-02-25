### Fitting VBGF Parameters for P. filametnosus using a JAGS based bayesian framework
proj_dir = getwd()
data_dir = file.path(proj_dir, "data")
src_dir = file.path(proj_dir, 'src')
results_dir = file.path(proj_dir, 'results')

### Installing principle dependencies
# install.packages("R2jags")
library('R2jags')
library('lattice')
library('coda')

### Loading in saved workspace
if(file.exists(file.path(results_dir, 'Bayes/Bayesian Workspace.RData'))){
  load(file.path(results_dir, 'Bayes/Bayesian Workspace.RData'))
}

### Loading in Datafile
mark_recapture_data = read.csv(file.path(data_dir, 'HO Mstr, temp (version 1).csv'), stringsAsFactors =  FALSE)
print(colnames(mark_recapture_data))
head(mark_recapture_data)

#### Cleaning Data
min_time_at_lib = 60

### Renaming data columns
colnames(mark_recapture_data) = c('tag_date', 'location', 'station', 'depth_f', 'species', 'previously_tagged', 'tag_id','fork_length_in', 'remarks', 'recapture_1_date', 'recapture_1_location', 'recapture_1_station', 'recapture_1_depth_f', 'recapture_1_fork_length_in', 'weight_1_lbs', 'days_1_free', 'growth_1_in', 'distance_1_miles','retagged_1',
                                  'recapture_2_date', 'recapture_2_location', 'recapture_2_station', 'recapture_2_depth_f', 'recapture_2_fork_length_in', 'weight_2_lbs', 'days_2_free', 'growth_2_in', 'distance_2_miles', 'retagged_2',
                                  'recapture_3_date', 'recapture_3_location', 'recapture_3_station', 'recapture_3_depth_f', 'recapture_3_fork_length_in', 'weight_3_lbs', 'days_3_free', 'growth_3_in', 'distance_3_miles', 'retagged_3',
                                  'recapture_4_date', 'recapture_4_location', 'recapture_4_station', 'recapture_4_depth_f', 'recapture_4_fork_length_in', 'weight_4_lbs', 'days_4_free', 'growth_4_in', 'distance_4_miles', 'x_retagged')


## How many total fish do we have in the data set?
n_tagged_fish = dim(mark_recapture_data)[1] # 4245!

### Subsetting out only Opakapaka with tag IDs - That is, fish that were marked
mark_recapture_data = mark_recapture_data[mark_recapture_data$species == '1' & mark_recapture_data$tag_id != '', ]
n_tagged_paka = dim(mark_recapture_data)[1] # This gets you  the previously published 4179 tagged paka number from Kobayashi, Okamoto, & Oishi . for some reason doesn't exclude fish marked 'died'

#### Adusting Data Classes
### Formatting Dates (Converting Characters to POSIXct)
mark_recapture_data$tag_date = as.POSIXct(mark_recapture_data$tag_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_1_date = as.POSIXct(mark_recapture_data$recapture_1_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_2_date = as.POSIXct(mark_recapture_data$recapture_2_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_3_date = as.POSIXct(mark_recapture_data$recapture_3_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_4_date = as.POSIXct(mark_recapture_data$recapture_4_date, format = "%Y-%m-%d")

### Formatting fork lengths
## We need to convert our forklength measurements from inches to cm and make sure everything is a numeric value.
## Note: There are a couple fork lengths that have ?, *, or have no lengths recorded. 
## I have no idea what these are but they're qualifiers and so I'm going to let them go to NA and get dropped from analysis
in_to_cm = 2.54
mark_recapture_data$fork_length_cm = as.numeric(mark_recapture_data$fork_length_in) * in_to_cm
mark_recapture_data$recapture_1_fork_length_cm = as.numeric(mark_recapture_data$recapture_1_fork_length_in) * in_to_cm
mark_recapture_data$recapture_2_fork_length_cm = as.numeric(mark_recapture_data$recapture_2_fork_length_in) * in_to_cm
mark_recapture_data$recapture_3_fork_length_cm = as.numeric(mark_recapture_data$recapture_3_fork_length_in) * in_to_cm
mark_recapture_data$recapture_4_fork_length_cm = as.numeric(mark_recapture_data$recapture_4_fork_length_in) * in_to_cm

### Removing fish that were not recaptured
mark_recapture_data = mark_recapture_data[!is.na(mark_recapture_data$recapture_1_fork_length_cm), ]
### This leaves us with 
dim(mark_recapture_data)[1] # individuals

#### Creating Model Data

### L: A matrix of fork lengths (cm). Rows are individuals, columns are capture events
L = cbind(mark_recapture_data$fork_length_cm, mark_recapture_data$recapture_1_fork_length_cm, mark_recapture_data$recapture_2_fork_length_cm, mark_recapture_data$recapture_3_fork_length_cm, mark_recapture_data$recapture_4_fork_length_cm)

### dt: A matrix of delta t (years) corrosponding to time between capture events. Column 1 is a dummy column that will be removed before data is wrapped in a list
dt = cbind(rep(9999, length(mark_recapture_data$recapture_1_date)), difftime(mark_recapture_data$recapture_1_date, mark_recapture_data$tag_date, unit = 'days') / 365,
           difftime(mark_recapture_data$recapture_2_date, mark_recapture_data$tag_date, unit = 'days') / 365,
           difftime(mark_recapture_data$recapture_3_date, mark_recapture_data$tag_date, unit = 'days') / 365,
           difftime(mark_recapture_data$recapture_4_date, mark_recapture_data$tag_date, unit = 'days') / 365)

### Ommitting data when the time at liberty is less than our defined minimum period
L[dt < min_time_at_lib / 365] = NA
dt[dt < min_time_at_lib / 365] = NA

### Removing any data from fish without a valid recapture event.
rm_ind = c()
## Loop through each individual
for(i in 1:nrow(L)){
  ## Flag individuals with less than 2 valid recaptures
  if(length(which(is.na(L[i, ]))) >= length(L[i, ]) - 1){
    rm_ind = c(rm_ind, i)
  } else if(length(which(is.na(dt[i, ]))) >= length(dt[i, ]) - 1){
    rm_ind = c(rm_ind, i)
  }
}
## Remove data flagged for removal
L = L[-rm_ind, ]; dt = dt[-rm_ind, ]


### Left justifying matricies dt and L
## Loop through each fish
for(i in 1:nrow(L)){
  ## If they have any NA values (total captures < 5)
  if(any(is.na(L[i, ]))){
    ## For as long as there is an NA value with a numeric right ajacent
    while(min(which(is.na(dt[i, ]))) < max(which(!is.na(dt[i, ])))){
      ## For matricies dt and L, create a new vector without the NA index, add a new NA to the rightmost side, and replace the current line 
      dt[i, ] = c(dt[i, 1:min(which(is.na(L[i, ])))-1], dt[i, (min(which(is.na(L[i, ])))+1):length(L[i, ])], NA)
      L[i, ] = c(L[i, 1:min(which(is.na(L[i, ])))-1], L[i, (min(which(is.na(L[i, ])))+1):length(L[i, ])], NA)
    }
  }
}

#### We should now have the same number of rows as tagdat ####
nrow(tagdat) == nrow(L)

#### dt still has a column on the right that is full of bullshit values
### Removing the first column of dt
dt = dt[ ,-1]

#### Getting values of n, a vector where values represent the number of valid recaptures for individuals
n = rep(0, nrow(dt))
for(i in 1:length(n)){
  n[i] = length(L[i, ]) - length(which(is.na(L[i, ])))
}


## Print number of rows in L. Should be 387
dim(L)

#### Wrapping all our data into a list for our model
data = list(n = n, L = L, dt = dt, N = length(n))

##### Model Specification #####
#### Defining inits - a list of initial variables
## Initial values for each chain are stored in 3 lists with elements corrosponding to all variables to be initialized
inits1 = list('Linf_mu' = 60,   'k_mu' = 0.30)
inits2 = list('Linf_mu' = 60*2, 'k_mu' = 0.3*2)
inits3 = list('Linf_mu' = 60/2, 'k_mu' = 0.3/2)
## Creating a single list that contains the lists corrosponding to each chain's initial values
inits = list(inits1, inits2, inits3)

#### Running our models ####

## Running a test model that omits intermediate recaptures.
## Running this test model should give similar paramter estimates to maximum likelihood method. 
# model_t = jags(data, inits, 
#                model.file = file.path(src_dir,'Bayes/ModelTest_Jags.txt'),
#                # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
#                parameters.to.save =  c('Linf_std', 'k_std', 'variance', 'Linf_mu', 'Linf_tau', 'Shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
#                DIC = T, 
#                n.chains = 3, n.iter = 10000, n.burnin = 100, n.thin = 50)
# 
# save.image(file.path(results_dir, 'Bayesian Workspace Img'))

## Model 1: Linf and K vary between individuals
model_1 = jags(data, inits, 
               model.file = file.path(src_dir,'Bayes/Model1_Jags.txt'),
               # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               parameters.to.save =  c('Linf_mu', 'Linf_std', 'Linf_tau', 'Shape', 'k_mu', 'k_std', 'k_tau', 'rate', 'tau', 'variance'),
               DIC = T, 
               n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50)
save.image(file.path(results_dir, 'Bayes/Bayesian Workspace.RData'))


## Model 2: Linf varies between individuals, K is fixed
model_2 = jags(data, inits, 
               model.file = file.path(src_dir,'Bayes/Model2_Jags.txt'),
               # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               parameters.to.save =  c('Linf_mu', 'Linf_std', 'Linf_tau', 'Shape', 'k_mu', 'k_std', 'k_tau', 'rate', 'tau', 'variance'),
               DIC = T, 
               n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50)
save.image(file.path(results_dir, 'Bayes/Bayesian Workspace.RData'))

## Model 3: K varies between individuals, Linf is fixed
model_3 = jags(data, inits, 
               model.file = file.path(src_dir,'Bayes/Model3_Jags.txt'),
               # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               parameters.to.save =  c('Linf_mu', 'Linf_std', 'Linf_tau', 'Shape', 'k_mu', 'k_std', 'k_tau', 'rate', 'tau', 'variance'),
               DIC = T, 
               n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50)
save.image(file.path(results_dir, 'Bayes/Bayesian Workspace.RData'))

## Model 4: Linf and K are fixed
model_4 = jags(data, inits, 
               model.file = file.path(src_dir,'Bayes/Model4_Jags.txt'),
               # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               parameters.to.save =  c('Linf_mu', 'Linf_std', 'Linf_tau', 'Shape', 'k_mu', 'k_std', 'k_tau', 'rate', 'tau', 'variance'),
               DIC = T, 
               n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50)
save.image(file.path(results_dir, 'Bayes/Bayesian Workspace.RData'))

growth_models = list('model_1' = model_1, 'model_2' = model_2, 'model_3' = model_3, 'model_4' = model_4)

##### Model Diagnostics
for(i in 1:length(growth_models)){
  # readline(paste('Press Enter to View Diagnostics for', names(growth_models)[i]))
  model = growth_models[[i]]
  
  ## Write out model summary table
  write.csv(model$BUGSoutput[10], file.path(results_dir, paste('Bayes/', names(growth_models)[i], ' Parameter Summaries.csv')))
}

for(i in 1:length(growth_models)){
  ### Summary Statistics and Parameter Estimates
  summary(model)
  
  ### Generating MCMC object for analysis
  model_mcmc = as.mcmc(model)
  
  ### xy plot
  xyplot(model_mcmc)
  
  ### Trace and Density Plots
  plot(model_mcmc)
  
  ### Autocorrelation plots
  autocorr.plot(model_mcmc)
  
  #### Other Diagnostics using CODA package
  gelman.plot(model_mcmc)
  geweke.diag(model_mcmc)
  geweke.plot(model_mcmc)
  raftery.diag(model_mcmc)
  heidel.diag(model_mcmc)
}


#### Extracting coefficients of variation for Linf and K paramters
cv_df = data.frame(stringsAsFactors = F)
for(i in 1:length(growth_models)){
  curr_mod = growth_models[[i]]
  cv_df = rbind(cv_df, data.frame('model' = names(growth_models)[i], 'Parameter' = 'Linf', 'cv' = curr_mod$BUGSoutput$summary["Linf_mu",'sd'] / curr_mod$BUGSoutput$summary["Linf_mu","mean"] * 100), 
                data.frame('model' = names(growth_models)[i], 'Parameter' = 'k', 'cv' = curr_mod$BUGSoutput$summary["k_mu",'sd'] / curr_mod$BUGSoutput$summary["k_mu","mean"] * 100))
}

library(forcats)
cv_df$`source of individual variability` = 'Other'
cv_df$`source of individual variability`[cv_df$model == 'model_1'] = 'Both'
cv_df$`source of individual variability`[cv_df$model == 'model_4'] = 'Neither'
cv_df$`source of individual variability`[cv_df$model == 'model_2' & cv_df$Parameter == 'Linf'] = 'Self'
cv_df$`source of individual variability`[cv_df$model == 'model_3' & cv_df$Parameter == 'k'] = 'Self'
cv_df$`source of individual variability` = as.factor(cv_df$`source of individual variability`)
cv_df$`source of individual variability` = fct_relevel(cv_df$`source of individual variability`, c('Both', 'Self', 'Other', 'Neither'))

library(ggplot2)
fig2 = ggplot(cv_df, aes(x = `source of individual variability`, y = cv, col = Parameter)) + geom_point() + geom_line(aes(group = Parameter)) + labs(x = 'Source of Individual Variability', y = 'Coefficient of Variation (Percent)', fill = 'Parameter') + theme(legend.justification = c(1, 1), legend.position = c(.15, .95)) + ylim(0,100)
pdf(file.path(results_dir, 'Bayes/Fig2 - Coefficients of Variation.pdf'))
  print(fig2)
dev.off()
                          