#### Utility Functions required for Comparing Age and Growth Increments from Bayesian and integrative data approaches for the deepwater snapper Pristipomoides filamentosus in the Hawaiian Islands

load_okomoto_data = function(){
  
  # A function for loading and cleaning the OTP dataset
  
  ## Load data file
  mark_recapture_data = read.csv(file.path(data_dir, 'HO Mstr, temp (version 1).csv'), stringsAsFactors =  FALSE)
  colnames(mark_recapture_data) = c('tag_date', 'location', 'station', 'depth_f', 'species', 'previously_tagged', 'tag_id','fork_length', 'remarks', 'recapture_1_date', 'recapture_1_location', 'recapture_1_station', 'recapture_1_depth_f', 'recapture_1_fork_length', 'weight_1_lbs', 'days_1_free', 'growth_1', 'distance_1_miles','retagged_1',
                                    'recapture_2_date', 'recapture_2_location', 'recapture_2_station', 'recapture_2_depth_f', 'recapture_2_fork_length', 'weight_2_lbs', 'days_2_free', 'growth_2', 'distance_2_miles', 'retagged_2',
                                    'recapture_3_date', 'recapture_3_location', 'recapture_3_station', 'recapture_3_depth_f', 'recapture_3_fork_length', 'weight_3_lbs', 'days_3_free', 'growth_3', 'distance_3_miles', 'retagged_3',
                                    'recapture_4_date', 'recapture_4_location', 'recapture_4_station', 'recapture_4_depth_f', 'recapture_4_fork_length', 'weight_4_lbs', 'days_4_free', 'growth_4', 'distance_4_miles', 'x_retagged')

  ## Subsetting out only fish with tag IDs (fish that were marked)
  mark_recapture_data = mark_recapture_data[mark_recapture_data$species == '1' & mark_recapture_data$tag_id != '', ]
  
  ## Format Dates 
  date_columns = c('tag_date', 'recapture_1_date', 'recapture_2_date', 'recapture_3_date', 'recapture_4_date')
  for (column in date_columns){
    mark_recapture_data[ , colnames(mark_recapture_data) == column] = as.POSIXct(mark_recapture_data[ , colnames(mark_recapture_data) == column], format = "%Y-%m-%d")
  }
 
  ## Transforming Fork Length from imperial to metric
  in_to_cm = 2.54
  fork_length_columns = c('fork_length', 'recapture_1_fork_length', 'recapture_2_fork_length', 'recapture_3_fork_length', 'recapture_4_fork_length')
  for (column in fork_length_columns){
    mark_recapture_data[ , colnames(mark_recapture_data) == column] = as.numeric(mark_recapture_data[ , colnames(mark_recapture_data) == column]) * in_to_cm
  }
  
  return(mark_recapture_data)
}


# Modifying our JOINT LIKELIHOOD function to accept more than one set of a given data type:
joint.logl.f <- function(param, npf, npA, tagdat = NULL, tagdat2 = NULL ,otodat = NULL, otodat2 = NULL, otodat3 = NULL, otodat4 = NULL, lfdat = NULL, lfdat2 = NULL, wt.oto=0, wt.oto2=0, wt.oto3=0, wt.oto4=0, wt.tag=0, wt.tag2=0,wt.lf=0, wt.lf2=0)
{

  neglogl.tag<- 0
  neglogl.tag2<- 0
  neglogl.oto<- 0
  neglogl.oto2<- 0
  neglogl.oto3<- 0
  neglogl.oto4<- 0
  neglogl.lf<- 0
  neglogl.lf2<- 0
  
  if(wt.tag>0)  neglogl.tag <- logl.ssnl.f(param,npf,npA,tagdat)
  if(wt.tag2>0)  neglogl.tag2 <- logl.ssnl.f(param,npf,npA,tagdat2)
  if(wt.oto>0)  neglogl.oto <- logl.oto.f(param,npf,npA,otodat)
  if(wt.oto2>0)  neglogl.oto2 <- logl.oto.f(param,npf,npA,otodat2)
  if(wt.oto3>0)  neglogl.oto3 <- logl.oto.f(param,npf,npA,otodat3)
  if(wt.oto4>0)  neglogl.oto4 <- logl.oto.f(param,npf,npA,otodat4)
  if(wt.lf>0) neglogl.lf <- logl.lf.f(param,npf,npA,lfdat)
  if(wt.lf2>0) neglogl.lf2 <- logl.lf.f(param,npf,npA,lfdat2)
  
  neglogl <- wt.tag*neglogl.tag + wt.tag2*neglogl.tag2 + wt.oto*neglogl.oto+ wt.oto2*neglogl.oto2 + wt.oto3*neglogl.oto3 + wt.oto4*neglogl.oto4 + wt.lf*neglogl.lf + wt.lf2*neglogl.lf2
  # print(c(neglogl.tag, neglogl.tag2, neglogl.oto, neglogl.oto2, neglogl.oto3, neglogl.oto4, neglogl.lf, neglogl.lf2, neglogl))
  print(neglogl)
  return(neglogl)
}

## Predicting length at recapture
predict_recapture_length = function(Lm, dt, linf = 65.95546, k = 0.2369113, a = 0){
  ## Get estimated length at recapture of a given individual using von Bertalanffy function as paramterized by Faben
  #return(linf * ((1 - exp(-k * (a + dt))) - (1 - exp(-k * a))))
  return(Lm + (linf - Lm) * (1 - exp(-k * dt)))
}

calculate_model_rmse = function(Lm, dt, linf, k, Lr_obs){
  return(sum((predict_recapture_length(Lm = Lm, dt = dt, linf = linf, k = k) - Lr_obs)^2) / length(Lr_obs))
}

std_error = function(x){
  #### Calculates standard error of set (x)
  sqrt(var(x)/length(x))
}

#### Determining how long it takes to reach a threshold % of Linf under each model
yrs_to_.9_linf = function(linf, k, a0 = 0, threshold = 0.90){
  t = log(1 - (linf * threshold/linf)) /  (-1 * k) + a0
  return(t)
}

#### Bootstrapping functions
### A function to bootstrap length-frequency data
lf_boot = function(pseudo_data){
  set.seed(Sys.time())
  ## For each month bin, replace with resampling the number of fish caught that month
  boot_lf_dat = NULL
  if(!exists('pseudo_data$curr_month_year')){
    pseudo_data$curr_month_year = pseudo_data$month_year
  }
  for(i in 1:length(unique(pseudo_data$curr_month_year))){
    monthly_pseudo_data = pseudo_data[pseudo_data$curr_month_year == unique(pseudo_data$curr_month_year)[i], ]
    monthly_n_fish = dim(monthly_pseudo_data)[1]
    boot_lf_dat = rbind(boot_lf_dat, monthly_pseudo_data[sample(x = 1:monthly_n_fish, size = monthly_n_fish, replace = TRUE), ])
  }
  return(boot_lf_dat)
}

## A function for bootstrap sampling with replacement
bootstrap_growth = function(boot_iterations = 10000, tagdat = NULL, tagdat2 = NULL, otodat = NULL, otodat2 = NULL, otodat3 = NULL, otodat4 = NULL, pseudolf = NULL, pseudolf2 = NULL, wt.oto = 1, wt.oto2 = 0, wt.oto3 = 0, wt.oto4 = 0, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0){
  boot_param_ests = NULL
  boot_param_ests = foreach(1:boot_iterations, .combine = rbind) %dopar%{
    boot_parms = NULL
    ## Resampling data
    if(!is.null(tagdat)){boot.tagdat = tagdat[sample(nrow(tagdat), size = nrow(tagdat), replace = TRUE), ]}else{boot.tagdat = NULL}
    if(!is.null(tagdat2)){boot.tagdat2 = tagdat[sample(nrow(tagdat2), size = nrow(tagdat2), replace = TRUE), ]}else{boot.tagdat2 = NULL}
    
    if(!is.null(otodat)){boot.otodat = otodat[sample(nrow(otodat), size = nrow(otodat), replace = TRUE), ]}else{boot.otodat = NULL}
    if(!is.null(otodat2)){boot.otodat2 = otodat2[sample(nrow(otodat2), size = nrow(otodat2), replace = TRUE), ]}else{boot.otodat2 = NULL}
    if(!is.null(otodat3)){boot.otodat3 = otodat3[sample(nrow(otodat3), size = nrow(otodat3), replace = TRUE), ]}else{boot.otodat3 = NULL}
    if(!is.null(otodat4)){boot.otodat4 = otodat4[sample(nrow(otodat4), size = nrow(otodat4), replace = TRUE), ]}else{boot.otodat4 = NULL}
    
    boot.lfdat = NULL
    if(!is.null(pseudolf)){
      while(class(boot.lfdat) != 'data.frame'){
        boot.lfdat  = try(length_freq_decomp(lf_boot(pseudolf)), silent = TRUE)
      }
    }
    
    boot.lfdat2 = NULL
    if(!is.null(pseudolf2)){
      while(class(boot.lfdat2) != 'data.frame'){
        boot.lfdat2  = length_freq_decomp(lf_boot(pseudolf))
      }
    }
    
    ## Refitting MLE Object
    boot.fit.vb = try(boot.fit.vb <- nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=boot.tagdat, tagdat2 = boot.tagdat2 ,otodat=boot.otodat, otodat2 = boot.otodat2, otodat3 = boot.otodat3, otodat4 = boot.otodat4, lfdat=boot.lfdat, lfdat2 = boot.lfdat2, wt.oto=wt.oto, wt.oto2=wt.oto2, wt.oto3=wt.oto3, wt.oto4=wt.oto4, wt.tag=wt.tag, wt.tag2=wt.tag2, wt.lf=wt.lf, wt.lf2=wt.lf2), silent = TRUE)
    ## Writing out parameters if model converged, otherwise writing out NAs
    if(is.list(boot.fit.vb)){
      return(boot.fit.vb$par)
    } else {
      return(rep(NA, 10))
    }
  }
  boot_param_ests = as.data.frame(boot_param_ests)
  colnames(boot_param_ests) = c('mu.L', 'sig.L',   'k',  'mu.A', 'sig.A', 'sig.sci', 'sig.f', 'a0', 'sig.oto', 'sig.lf')
  return(boot_param_ests)
}

## A function for producing summary statistics from the bootstrap parameter distributions
calc_boot_stats = function(boot_param_ests){
  boot_param_ests = boot_param_ests[!is.na(boot_param_ests[,1]), ]
  two.five.percent = dim(boot_param_ests)[1] * 0.025
  ninetyseven.five.percent = dim(boot_param_ests)[1] * 0.975
  boot_stats = data.frame(matrix(NA, nrow = dim(boot_param_ests)[2], ncol = 3), stringsAsFactors = FALSE)
  colnames(boot_stats) = c('Median', '2.5%', '97.5%')
  rownames(boot_stats) = colnames(boot_param_ests)
  for(i in 1:dim(boot_param_ests)[2]){
    boot_stats[i, 1] = median(boot_param_ests[ ,i], na.rm = TRUE)
    if(floor(two.five.percent) != two.five.percent){
      boot_stats[i, 2] = (sort(boot_param_ests[ ,i])[floor(two.five.percent)] * abs(two.five.percent - floor(two.five.percent))) + (sort(boot_param_ests[ ,i])[ceiling(two.five.percent)] * abs(two.five.percent - ceiling(two.five.percent)))
    } else {
      boot_stats[i, 2] = sort(boot_param_ests[ ,i])[floor(two.five.percent)]
    }
    if(floor(ninetyseven.five.percent) != ninetyseven.five.percent){
      boot_stats[i, 3] = (sort(boot_param_ests[ ,i])[floor(ninetyseven.five.percent)] * abs(ninetyseven.five.percent - floor(ninetyseven.five.percent))) + (sort(boot_param_ests[ ,i])[ceiling(ninetyseven.five.percent)] * abs(ninetyseven.five.percent - ceiling(ninetyseven.five.percent)))
    } else {
      boot_stats[i, 3] = sort(boot_param_ests[ ,i])[floor(ninetyseven.five.percent)]
    }
  }
  return(t(boot_stats))
}

## A function to bring all bootstrapped function results together
bootstrap_growth_params = function(boot_iterations = 10000, filename = NULL, tagdat = NULL, tagdat2 = NULL, otodat = NULL, otodat2 = NULL, otodat3 = NULL, otodat4 = NULL, pseudolf = NULL, pseudolf2 = NULL, wt.oto = 1, wt.oto2 = 0, wt.oto3 = 0, wt.oto4 = 0, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0){
  if(boot_iterations < 40){
    print('Error: boot_iterations must be > 40')
    return()
  }
  ## First run initial bootstrapping pass
  booted_param_ests = bootstrap_growth(boot_iterations, tagdat, tagdat2, otodat, otodat2, otodat3, otodat4, pseudolf, pseudolf2, wt.oto, wt.oto2, wt.oto3, wt.oto4, wt.tag, wt.tag2, wt.lf, wt.lf2)
  
  ## Then rerun to fill in iterations that failed to converge
  while(dim(booted_param_ests[!is.na(booted_param_ests$mu.L), ])[1] < boot_iterations){
    booted_param_ests = rbind(booted_param_ests, bootstrap_growth(boot_iterations = length(which(is.na(booted_param_ests$mu.L))), tagdat, tagdat2, otodat, otodat2, otodat3, otodat4, pseudolf, pseudolf2, wt.oto, wt.oto2, wt.oto3, wt.oto4, wt.tag, wt.tag2, wt.lf, wt.lf2))
  }
  
  ## What percentage failed to converge?
  convergence_failure_rate = (length(which(is.na(booted_param_ests$mu.L))) / dim(booted_param_ests)[1]) * 100
  print(paste('Convergence failure rate =', round(convergence_failure_rate, digits = 2), '%'))
  
  ## Getting summary stats and writing them out
  boot_stats = calc_boot_stats(booted_param_ests)
  if(!is.null(filename)){
    write.csv(boot_stats, filename)
    write.csv(booted_param_ests,filename)
  }
  
  ## Writing out results
  results = list()
  results$raw_boot_data = booted_param_ests
  results$boot_stats = boot_stats
  results$convergence_failure_rate = convergence_failure_rate
  return(results)
}



len = length



### Age at recruitment to juvenile fishing grounds. 
## Peak spawning occurs in July according to Luers, Demartini, Humphreys 2017
## Difference between first sampling trip (Oct) and peak spawning (July) = 3 months




## Creating a function that we can use later for bootstrapping
length_freq_decomp = function(pseudo_data, plot = FALSE, fixed_modes = FALSE){
  # Function for decompositing length frequency data
  
  ##  starting means estimated from data
  start_means = rbind(c(1, 17), c(11, 18), c(13, 20), c(14, 20), c(15, NA), c(15, NA), c(16, NA), c(17, NA), c(16, NA), c(10, 16), c(13, 16), c(13, 18), c(15, 17))
  constrain_means = rbind(c(11, 17),c(11, 18), c(14, 19), c(13, 19), c(14, NA), c(14, NA), c(14, NA), c(17, NA), c(15, NA), c(10, 16.5), c(13, 17), c(13, 17), c(15, 19))
  
  ### Age at recruitment to juvenile fishing grounds. 
  ## Peak spawning occurs in July according to Luers, Demartini, Humphreys 2017
  ## Difference between first sampling trip (Oct) and peak spawning (July) = 3 months
  age_at_recruitment = 3/12
  
  lfdat = data.frame(stringsAsFactors = FALSE)
  for(i in 1:length(unique(pseudo_data$date))){
    curr_month_year = as.character(unique(pseudo_data$date)[i])
    month_data = pseudo_data[pseudo_data$date == curr_month_year, ]
    ## Determing the mean age of fish once they've recruited, assuming that they were born during peak spawning.
    mode.age = as.numeric(difftime(unique(month_data$date), as.POSIXct('1989-07-01'), "days")) / 365 # in years
    mode.age = c(mode.age, mode.age + 1) # Assumption is that if two cohorts are present, the second is one year older than the first
    # When we get to the second year of data, YOY for first year becomes 1+ year old, new cohort recruits. Because we based age on difftime for first cohort, we need to remove 1 year from all ages
    if(min(mode.age) > (1 + age_at_recruitment)){ 
      mode.age = mode.age - 1
    }
    
    ## We need to decompose distributions for each cohort present into age of cohort, mean length of cohort, se of length of cohort, and number of fish falling into each cohort
    ## During Oct - Feb, two age cohorts present, we need to decompose two normal distributions from data (Moffitt and Parrish 1996)
    
    if (months(unique(month_data$date)) %in% c('October', 'November', 'January', 'February')) {
      decomp = NULL
      while(class(decomp) != 'mixEM'){
        k = 2 # Number of cohorts
        if(fixed_modes == TRUE){
          decomp = try(normalmixEM(month_data$length, mu = start_means[i, ], k = k, arbvar = TRUE, mean.constr = constrain_means[i, ]), silent = TRUE) # abvar = FALSE would both cohorts to have same sigma. See justification for this in Laslett et. al 2004
        } else {
          decomp = try(normalmixEM(month_data$length, mu = start_means[i, ], k = k, arbvar = TRUE), silent = TRUE) # abvar = FALSE would both cohorts to have same sigma. See justification for this in Laslett et. al 2004
        }
      }
      mode.len = decomp$mu[order(decomp$mu)] # Sometimes things pop out in a weird order. We assume that the smaller size class is the younger cohort
      est.n = c(decomp$lambda * dim(month_data)[1])[order(decomp$mu)]
      mode.se = (decomp$sigma / sqrt(est.n))[order(decomp$mu)]
      lfdat = rbind(lfdat, cbind(mode.age, mode.len, mode.se, est.n, curr_month_year))
    } else {
      k = 1 # Number of cohorts
      est.n = length(month_data$length) # Number of fish in cohort
      
      ## For times when only a single cohort exists, we need to figure out the relative age of that cohort
      if (format(month_data$date, "%m")[1] < 7) { # if month is less than october (month that new recruits show up) 
        mode.age = mode.age[1] # Go with the younger year class because fish are less than 1 year old
      } else {
        mode.age = mode.age[2] # Go with the older year class because fish are older than 1 (YOY have not recruited yet)
      }
      
      mode.len = mean(month_data$length)
      mode.se = sd(month_data$length)  / sqrt(est.n)
      
      if(fixed_modes == TRUE){
        mode.len = constrain_means[i, ][which(!is.na(constrain_means[i, ]))]
        mode.se = sqrt(sum((mode.len - month_data$length)^2)/(est.n - 1)) / est.n
      }
      
      lfdat = rbind(lfdat, cbind(mode.age, mode.len, mode.se, est.n, curr_month_year))
    }
  }
  
  ## We may need to reclass our data depending on if we wrote curr_month_year into lfdat or not (this turns stuff into factors)
  if(is.factor(lfdat$mode.age)){
    lfdat$mode.age = as.numeric(levels(lfdat$mode.age)[lfdat$mode.age])
    lfdat$mode.len = as.numeric(levels(lfdat$mode.len)[lfdat$mode.len])
    lfdat$mode.se  = as.numeric(levels(lfdat$mode.se)[lfdat$mode.se])
  } else {
    lfdat$mode.age = as.numeric(lfdat$mode.age)
    lfdat$mode.len = as.numeric(lfdat$mode.len)
    lfdat$mode.se  = as.numeric(lfdat$mode.se)
  }
  
  ## Sorting lfdat by the mean age. This makes it easier to visually inspect that fish are getting larger as they get older. This can get messed up during bimodal composition. 
  lfdat = lfdat[order(lfdat$mode.age), ]
  if(plot){
    plot(x = lfdat$mode.age, y = lfdat$mode.len, pch = 19)
  }
  return(lfdat)
}


calculate_model_rmse = function(Lm, dt, linf, k, Lr_obs){
  return(sqrt(sum((predict_recapture_length(Lm = Lm, dt = dt, linf = linf, k = k) - Lr_obs)^2) / length(Lr_obs)))
}

evaluate_models = function(cross_validation_iterations = 10000){
  mod_eval_results = data.frame(stringsAsFactors = FALSE)
  model_na_results = rep(0, 7)
  lit_vbgf = lit_vbgc_params[lit_vbgc_params$region %in% c('Hawaii - MHI & NWHI', 'Hawaii - MHI', 'Hawaii - NWHI'), ]
  lit_vbgf_for_train = lit_vbgf[!(lit_vbgf$author %in% paste('Bayesian Model',1:4)), ]
  bayes_models = lit_vbgf[(lit_vbgf$author %in% paste('Bayesian Model',1:4)), ]
  
  mod_eval_results = foreach(i = 1:cross_validation_iterations, .combine = rbind) %dopar% {
    ## Train Test Split
    train_index = sample(1:dim(tagdat)[1], size = n_train, replace = FALSE)
    tagdat_train = tagdat[train_index, ]
    # train2_index = sample(1:dim(tagdat2)[1], size = n_train2, replace = FALSE)
    # tagdat2_train = tagdat2[train2_index, ]
    tagdat_test = tagdat[-train_index, ]#rbind(tagdat[-train_index, ], tagdat2[-train2_index, ])
    
    score = rep(NA, 7)
    ### Setting intial params for all data
     #        mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
     p0 <- c(  70,     1,  .10,     1,   .10,      1,     0,    .3,      1,      1)
     lb <- c(  40,   0.01, .05,  0.01,   .05,   0.01,     0,  -.4,   0.01,   0.01)
     ub <- c( 110,  15.0,  .50,   1.5,   .50,   15.0,     0,   -.2,     15,     15)

    var5 = NULL
    var5 = try(nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=tagdat_train, tagdat2 = tagdat2_train,  wt.oto=0,wt.tag=1, wt.tag2 = 0, wt.lf=0)$par, silent = TRUE)
    if(is.numeric(var5[1])){
      score[1] = calculate_model_rmse(Lm = tagdat_test[ ,1], dt = tagdat_test[ ,4], linf = var5[1], k = var5[3], Lr_obs = tagdat_test[ ,2])
    } else {
      model_na_results[1] = model_na_results[1] + 1
    }
    
     ub <- c( 110,  15.0,  10,   1.5,   .50,   15.0,     0,   15,     15,     15)
    var6 = NULL
    var6 = try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = tagdat2_train, otodat=otodat, lfdat=lfdat, wt.oto=1/dim(otodat)[1], wt.tag=1/dim(tagdat_train)[1], wt.tag2=0, wt.lf=1/length(lfdat$curr_month_year))$par, silent = TRUE)
    if(is.numeric(var6[1])){
      score[2] = calculate_model_rmse(Lm = tagdat_test[ ,1], dt = tagdat_test[ ,4], linf = var6[1], k = var6[3], Lr_obs = tagdat_test[ ,2])
    } else {
      model_na_results[2] = model_na_results[2] + 1
    }
    
    var7 = NULL
    var7 = try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = tagdat2_train, otodat=otodat, lfdat=lfdat, wt.oto=1, wt.tag=1, wt.tag2 = 0, wt.lf=1)$par, silent = TRUE)
    if(is.numeric(var7[1])){
      score[3] = calculate_model_rmse(Lm = tagdat_test[ ,1], dt = tagdat_test[ ,4], linf = var7[1], k = var7[3], Lr_obs = tagdat_test[ ,2])
    } else {
      model_na_results[3] = model_na_results[3] + 1
    }
    
    var8 = NULL
    var8 =  try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = tagdat2_train, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 1/dim(otodat[otodat$source == 'ralston and miyamoto', ])[1], wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat_train)[1], wt.tag2 = 0, wt.lf = 1/length(lfdat$curr_month_year), wt.lf2 = 0)$par, silent = TRUE)
    if(is.numeric(var8[1])){
      score[4] = calculate_model_rmse(Lm = tagdat_test[ ,1], dt = tagdat_test[ ,4], linf = var8[1], k = var8[3], Lr_obs = tagdat_test[ ,2])
    } else {
      model_na_results[4] = model_na_results[4] + 1
    }
    
    var9 = NULL
    var9 =  try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = tagdat2_train, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1, wt.oto2= 1, wt.oto3=1, wt.oto4=1, wt.tag = 1,wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)$par, silent = TRUE)
    if(is.numeric(var9[1])){
      score[5] = calculate_model_rmse(Lm = tagdat_test[ ,1], dt = tagdat_test[ ,4], linf = var9[1], k = var9[3], Lr_obs = tagdat_test[ ,2])
    } else {
      model_na_results[5] = model_na_results[5] + 1
    }
    
    var10 = NULL
    var10 =  try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = tagdat2_train, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 0, wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat_train)[1], wt.tag2=0, wt.lf = 1/length(lfdat$curr_month_year), wt.lf2 = 0)$par, silent = TRUE)
    if(is.numeric(var10[1])){
      score[6] = calculate_model_rmse(Lm = tagdat_test[ ,1], dt = tagdat_test[ ,4], linf = var10[1], k = var10[3], Lr_obs = tagdat_test[ ,2])
    } else {
      model_na_results[6] = model_na_results[6] + 1
    }
    
    var11 = NULL
    var11 = try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = tagdat2_train, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat = lfdat, lfdat2 = NULL, wt.oto = 1, wt.oto2= 0, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)$par, silent = TRUE)
    if(is.numeric(var11[1])){
      score[7] = calculate_model_rmse(Lm = tagdat_test[ ,1], dt = tagdat_test[ ,4], linf = var11[1], k = var11[3], Lr_obs = tagdat_test[ ,2])
    } else {
      model_na_results[7] = model_na_results[7] + 1
    }
    
    ## Now getting fits from literature data  
    ## lit_models
    lit_var_scores = rep(0, length(lit_vbgf_for_train$author))
    for(i in 1:length(lit_vbgf_for_train$author)){
      lit_var_scores[i] = calculate_model_rmse(Lm = tagdat_test[ ,1], dt = tagdat_test[ ,4], linf = lit_vbgf_for_train$linf[i], k = lit_vbgf_for_train$k[i], Lr_obs = tagdat_test[ ,2])
    }
    
    ## Bayes_models
    bayes_models = lit_vbgc_params[21:24, ]
    bayes_var_scores = rep(0, length(bayes_models$author))
    for(i in 1:length(bayes_models$author)){
      bayes_var_scores[i] = calculate_model_rmse(Lm = tagdat_test[ ,1], dt = tagdat_test[ ,4], linf = bayes_models$linf[i], k = bayes_models$k[i], Lr_obs = tagdat_test[ ,2])
    }
    
    ## Best Overall Model
    return(c(score, lit_var_scores, bayes_var_scores))
  }
  colnames(mod_eval_results) = c(paste('model', 5:11), lit_vbgf_for_train$author, bayes_models$author)
  return(invisible(as.data.frame(mod_eval_results)))
}

# 
# gstats.ssnl.f<- function(tagdat,param.g,param.A,param.sig){   
#   qt1<- 0.001
#   qt2<- 15
#   
#   l<- matrix(rep(qt1,nrow(tagdat)),ncol=1)
#   u<- matrix(rep(qt2,nrow(tagdat)),ncol=1)
#   
#   aa<- gridsch.f(logint.ssnl.f,l,u,tagdat,param.g,param.A,param.sig)
#   brk<-0
#   for(i in (1:10))
#   {ms<- quadmax.f(aa)
#   
#   # we need to make sure that the mng values have stabilised
#   # if it has we break out
#   if(i>3)
#   {stdiff<- abs((ms[,1]-mngold)/ms[,2])
#   msd<- max(stdiff)
#   if(msd<0.2 & !is.nan(msd))
#   {   brk<- 1
#   break
#   }
#   }
#   
#   # if not, we refine the estimates
#   a2<- pmax(ms[,1],0.9*aa[,1]+0.1*aa[,5])
#   a2<- pmin(a2,0.1*aa[,1]+0.9*aa[,5])
#   a1<- pmin(a2,aa[,3])
#   a2<- a2+aa[,3]-a1
#   a1<- pmin(a1,0.1*aa[,1]+0.9*a2)
#   f12<- logint.ssnl.f(cbind(a1,a2),tagdat,param.g,param.A,param.sig)
#   tf<- f12[,1]<f12[,2]
#   ### ADDED LINE BY SS 27 May 2020
#   tf <- tf[!is.na(tf)]
#   ### 
#   atemp<- cbind(a1,f12[,1],a2,f12[,2],aa[,5:6])
#   aa<- cbind(aa[,1:2],atemp[,1:4])
#   aa[tf,]<-atemp[tf,]
#   mngold<- ms[,1]
#   }
#   if(brk==0) {assign("nits.gstats",1) #,where=1)
#     print("max # iterations reached")}
#   ms<- quadmax.f(aa)
#   return(ms)
# }
# 
# 
# log.lognormA.f<- function(a,param.A)
# {   # the log density of the tagging age A 
#   # A log-normal 
#   alpha.A<- param.A[1]
#   beta.A<-  param.A[2]
#   
#   # a is nf x n matrix
#   lna<- a
#   
#   ### THIS LINE Changed BY SS 27 MAY 2020
#   tf<- a>0 & !is.na(a)
#   
#   lna[tf]<- log(a[tf])
#   lna[!tf]<- -300
#   
#   logpa<- (lna-alpha.A)/beta.A
#   logpa<- -0.5*log(2*pi)-log(beta.A)-lna-0.5*logpa^2
#   return(logpa)
# }
# 
# logdensA.f<- log.lognormA.f
# densA.f<-lognormA.f
# 


#### THIS FUNCTION WORKS FOR MODELs 5 and 11
##       mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
p0 <- c(  70,   3, .10,     1,   .10,      1,     0,       1,       1,      1)

lb <- c(  60,   2.1,   .01,    .5,   .05,    0.1,    0,   -2.5,     0.5,     0)
ub <- c(  90,   5,    .4,     2,   0.8,      4,    0,      0,       8,     5)

generate_sensitive_data = function(linf, k){
  ## A wrapper function for generating synthetic data form a uniform distribution for Maximum Likelihood analysis
  tagdat$bin = floor(tagdat$L1 / 5) * 5
  
  sensitive_data = data.frame()
  
  length_bins = seq(15, max(tagdat$bin), by = 5)
  for (bin in length_bins){
    obs = length(tagdat$L1[tagdat$bin == bin])
    bin_sd = sd(tagdat$L1[tagdat$bin == bin])
    mean_dt = mean(tagdat$dt[tagdat$bin == bin])
    
    if(is.na(bin_sd) & bin < 25){
      bin_sd = sd(tagdat$L1[tagdat$bin == 20])
    }
    
    if(is.na(bin_sd) & bin > 50){
      bin_sd = sd(tagdat$L1[tagdat$bin == 50])
    }
    
    L1 = rnorm(200 - obs, mean = bin + 2.5, sd = 0.001)
    dt = abs(rnorm(200 - obs, mean = mean(tagdat$dt), sd = sd(tagdat$dt)))
    L2 = predict_recapture_length(Lm = L1, dt = dt, linf = linf, k = k)
    sensitive_data = rbind(sensitive_data, data.frame('L1' = L1, 'L2' = L2, 'X.' = 0, 'dt' = dt, 'L2measurer' = 0))
  }
  
  sensitive_data = rbind(sensitive_data, tagdat[ ,-dim(tagdat)[2]])
  
  bad_ind = which(sensitive_data$L1 >= sensitive_data$L2)
  x = sensitive_data$L1[bad_ind] 
  sensitive_data$L1[bad_ind] = sensitive_data$L2[bad_ind]
  sensitive_data$L2[bad_ind] = x
  
  # sensitive_data$L2[sensitive_data$L2 <= sensitive_data$L1] = sensitive_data$L1[sensitive_data$L2 < sensitive_data$L1] + 0.001
  
  return(sensitive_data)
}
