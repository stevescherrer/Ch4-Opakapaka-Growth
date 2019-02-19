##### Fitting VBGF for Pristipomoides filamentosus with Mark Recapture Data
#### Written by: Stephen Scherrer with some code modified from Erik Franklin (2017)
#### Written: Feb - March 2018
#### Contact: scherrer@hawaii.edu
#### All Wrongs Preserved

#### Laslett et al 2004 implementation of Fabens Method using tagging data as well as length frequency and direct ageing data

##### Workspace Setup #####
## Clearing workspace
rm(list = ls())
print('Opakapaka Growth Analysis')
initial_run_time = Sys.time()
print(initial_run_time)

## Setting a Script Timer
script_timer = proc.time()

## Declaring Directory Path
proj_dir = getwd()
if(!"Okamoto_Mark_Recapture" %in% strsplit(proj_dir, '/')[[1]]){
  proj_dir = file.path(getwd(), "Okamoto_Mark_Recapture")
}
data_dir = file.path(proj_dir, "data")
src_dir = file.path(proj_dir, "src")
results_dir = file.path(proj_dir, "results")

## Creating a run specific results folder
run_results_dir = file.path(results_dir, paste('run', initial_run_time))
dir.create(run_results_dir)

print(paste('proj_dir:', proj_dir))
print(paste('src_dir:', src_dir))
print(paste('run_results_dir:', run_results_dir))

## Installing Principle Dependencies
print('Installing principle dependencies')
 library('notifyR') # #send_push()
library('doParallel')
# library('beepr')
library('mixtools')

## Sourcing R Scripts provided by Eveson/Laslett
print('Sourcing files')
source(file.path(src_dir, "Laslett Functions/joint_lkhd.r"))
source(file.path(src_dir, "Laslett Functions/growth_functions.r"))
source(file.path(src_dir, "Laslett Functions/tag_lkhd.r"))


## Reading in literature parameter values
print('Loading Data')
lit_vbgc_params = read.csv(file.path(data_dir, "Parameter Estimates.csv"), stringsAsFactors = FALSE)
lit_vbgc_params = lit_vbgc_params[!is.na(lit_vbgc_params$Linf), ]
colnames(lit_vbgc_params) = c('author', 'n', 'linf', 'k', 't0', 'region', 'method')
lit_vbgc_params = lit_vbgc_params[c(1:20, 22:25), ]

## Assigning cores for parallel processing
registerDoParallel(cores = detectCores()-1)

##### Defining Utility Functions #####  
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
  # print(param)
  # print(c(neglogl.tag, neglogl.tag2, neglogl.oto, neglogl.oto2, neglogl.oto3, neglogl.oto4, neglogl.lf, neglogl.lf2, neglogl))
  return(neglogl)
}

## Predicting length at recapture
predict_recapture_length = function(Lm, dt, linf = 65.95546, k = 0.2369113, a = 0){
  ## Get estimated length at recapture of a given individual using von Bertalanffy function as paramterized by Faben
  #return(linf * ((1 - exp(-k * (a + dt))) - (1 - exp(-k * a))))
  return(Lm + (linf - Lm) * (1 - exp(-k * dt)))
}

calculate_predictive_variance = function(Lm, dt, linf, k, Lr_obs){
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
    write.csv(boot_stats, file.path(run_results_dir, paste(filename, '.csv', sep = "")))
    write.csv(booted_param_ests, file.path(run_results_dir, paste(filename, '_raw.csv', sep = "")))
  }
  
  ## Writing out results
  results = list()
  results$raw_boot_data = booted_param_ests
  results$boot_stats = boot_stats
  results$convergence_failure_rate = convergence_failure_rate
  return(results)
}


##### Loading and Cleaning Data Files #####
#### Mark Recapture Data ####
mark_recapture_data = read.csv(file.path(data_dir, 'HO Mstr, temp (version 1).csv'), stringsAsFactors =  FALSE)

### Renaming data columns
colnames(mark_recapture_data) = c('tag_date', 'location', 'station', 'depth_f', 'species', 'previously_tagged', 'tag_id','fork_length_in', 'remarks', 'recapture_1_date', 'recapture_1_location', 'recapture_1_station', 'recapture_1_depth_f', 'recapture_1_fork_length_in', 'weight_1_lbs', 'days_1_free', 'growth_1_in', 'distance_1_miles','retagged_1',
                                  'recapture_2_date', 'recapture_2_location', 'recapture_2_station', 'recapture_2_depth_f', 'recapture_2_fork_length_in', 'weight_2_lbs', 'days_2_free', 'growth_2_in', 'distance_2_miles', 'retagged_2',
                                  'recapture_3_date', 'recapture_3_location', 'recapture_3_station', 'recapture_3_depth_f', 'recapture_3_fork_length_in', 'weight_3_lbs', 'days_3_free', 'growth_3_in', 'distance_3_miles', 'retagged_3',
                                  'recapture_4_date', 'recapture_4_location', 'recapture_4_station', 'recapture_4_depth_f', 'recapture_4_fork_length_in', 'weight_4_lbs', 'days_4_free', 'growth_4_in', 'distance_4_miles', 'x_retagged')


## How many total fish do we have in the data set?
dim(mark_recapture_data)[1] # 4245!

### Subsetting out only Opakapaka with tag IDs - That is, fish that were marked
mark_recapture_data = mark_recapture_data[mark_recapture_data$species == '1' & mark_recapture_data$tag_id != '', ]
  dim(mark_recapture_data)[1] # This gets you to the previously published 4179 tagged paka number from Kobayashi, Okamoto, & Oishi . for some reason doesn't exclude fish marked 'died'

#### Adusting Data Classes
### Formatting Dates (Converting Characters to POSIXct)
mark_recapture_data$tag_date = as.POSIXct(mark_recapture_data$tag_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_1_date = as.POSIXct(mark_recapture_data$recapture_1_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_2_date = as.POSIXct(mark_recapture_data$recapture_2_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_3_date = as.POSIXct(mark_recapture_data$recapture_3_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_4_date = as.POSIXct(mark_recapture_data$recapture_4_date, format = "%Y-%m-%d")

### Formatting fork lengths 
## Note: There are a couple fork lengths that have ?, *, or have no lengths recorded. 
## I have no idea what these are but they're qualifiers and so I'm going to let them go to NA and get dropped from analysis
in_to_cm = 2.54
mark_recapture_data$fork_length_cm = as.numeric(mark_recapture_data$fork_length_in) * in_to_cm
mark_recapture_data$recapture_1_fork_length_cm = as.numeric(mark_recapture_data$recapture_1_fork_length_in) * in_to_cm
mark_recapture_data$recapture_2_fork_length_cm = as.numeric(mark_recapture_data$recapture_2_fork_length_in) * in_to_cm
mark_recapture_data$recapture_3_fork_length_cm = as.numeric(mark_recapture_data$recapture_3_fork_length_in) * in_to_cm
mark_recapture_data$recapture_4_fork_length_cm = as.numeric(mark_recapture_data$recapture_4_fork_length_in) * in_to_cm


#### Now we want to format a table with lm (length at marking), lr (length at recapture), and dt (elapsed time)
### Note: If a fish was recaptured multiple times, there is a single entry for that individual corrosponding to the length at initial marking and the length at final recapture
paka_growth = data.frame(stringsAsFactors = FALSE)
for(i in 1:length(mark_recapture_data$tag_id)){
  if(!is.na(mark_recapture_data$recapture_4_fork_length_cm[i])){
    paka_growth = rbind(paka_growth, data.frame('tag_id' = mark_recapture_data$tag_id[i], 'Lm' = mark_recapture_data$fork_length_cm[i], 'Lr' = mark_recapture_data$recapture_4_fork_length_cm[i], 'tm' = mark_recapture_data$tag_date[i], 'tr' = mark_recapture_data$recapture_4_date[i], 'n_recaptures' = 4, stringsAsFactors = FALSE))
  }else if(!is.na(mark_recapture_data$recapture_3_fork_length_cm[i])){
    paka_growth = rbind(paka_growth, data.frame('tag_id' = mark_recapture_data$tag_id[i], 'Lm' = mark_recapture_data$fork_length_cm[i], 'Lr' = mark_recapture_data$recapture_3_fork_length_cm[i], 'tm' = mark_recapture_data$tag_date[i], 'tr' = mark_recapture_data$recapture_3_date[i], 'n_recaptures' = 3, stringsAsFactors = FALSE))
  }else if(!is.na(mark_recapture_data$recapture_2_fork_length_cm[i])){
    paka_growth = rbind(paka_growth, data.frame('tag_id' = mark_recapture_data$tag_id[i], 'Lm' = mark_recapture_data$fork_length_cm[i], 'Lr' = mark_recapture_data$recapture_2_fork_length_cm[i], 'tm' = mark_recapture_data$tag_date[i], 'tr' = mark_recapture_data$recapture_2_date[i], 'n_recaptures' = 2, stringsAsFactors = FALSE))
  }else if(!is.na(mark_recapture_data$recapture_1_fork_length_cm[i])){
    paka_growth = rbind(paka_growth, data.frame('tag_id' = mark_recapture_data$tag_id[i], 'Lm' = mark_recapture_data$fork_length_cm[i], 'Lr' = mark_recapture_data$recapture_1_fork_length_cm[i], 'tm' = mark_recapture_data$tag_date[i], 'tr' = mark_recapture_data$recapture_1_date[i], 'n_recaptures' = 1, stringsAsFactors = FALSE))
  }
}
paka_growth$dt = abs(difftime(paka_growth$tm, paka_growth$tr, units = "days"))    ## Converting dt from days to years
paka_growth$dt = as.numeric(paka_growth$dt) / 365 # Converting to years
### Constructing derived variable dl (change in growth)
paka_growth$dL = paka_growth$Lr - paka_growth$Lm
### Removing any fish that have a NA value for dL or dt (There is a single fish which had no tagging length and 7 fish with no recapture dates)
length(which(is.na(paka_growth$dL))) # 1
length(which(is.na(paka_growth$dt))) # 7
paka_growth = paka_growth[!is.na(paka_growth$dL) & !is.na(paka_growth$dt), ]

#### Creating a subset data frame that removes recording errors in length and time
# paka_growth = subset(paka_growth, dL > 0)
length(which(paka_growth$dt <= 60/365)) #46
paka_growth = subset(paka_growth, dt >= 60/365)
tagdat = as.matrix(data.frame('L1' = paka_growth$Lm, "L2" = paka_growth$Lr, " " = rep(0, length(paka_growth$Lr)), "dt" = paka_growth$dt, "L2measurer" = rep(0, length(paka_growth$Lr))))

#### Splitting data into training and validation sets
n_train    = round(dim(tagdat)[1] * (2/3))
n_validate = dim(tagdat)[1] - n_train
train_index = sample(1:dim(tagdat)[1], size = n_train, replace = FALSE)
tagdat_train = tagdat[train_index, ]
tagdat_validate = tagdat[-train_index, ]

#### Otolith Data (Ralston and Miyamoto 1983, DeMartini et al. 1994, Andrews et al. 2012) ####
otodat = read.csv(file.path(data_dir, "RalstonMiyamotoandDemartiniAndrews.csv"))
colnames(otodat) = c("age", "len", "source")

#### Length Frequency Data #### 
### Extrapolated from Moffitt and Parrish 1996 using earlier version of manuscript (1994)
### Monthly fish length counts were extrapolated from histograms in paper first by adjusting rotation so histograms were 'level', then fitting equally spaced bars for each n across y axis and comparing bar heights

### I ended up with one more record than they do but pretty close! (1048 vs. 1047).  number next to each month is the total number of fish for that month estimated from the histograms. * means this number was estimated a second time and matched
oct_1989 = data.frame('date' = as.POSIXct('1989-10-01'), 'val' = c(0, 0, 1, 2, 1, 7, 1, 2, 1, 8, 7, 18, 5, 3, 0, 0, 0, 0, 0), len = 6:24) # 56 *
nov_1989 = data.frame('date' = as.POSIXct('1989-11-01'), 'val' = c(0, 0, 1, 4, 11, 7, 6, 4, 4, 3, 11, 12, 5, 1, 0, 0, 0, 0, 0), len = 6:24) # 69 *
jan_1990 = data.frame('date' = as.POSIXct('1990-01-01'), 'val' = c(0, 0, 0, 1, 6, 10, 12, 13, 20, 8, 1, 5, 2, 6, 3, 5, 1, 0, 1), len = 6:24) # 94 *
feb_1990 = data.frame('date' = as.POSIXct('1990-02-01'), 'val' = c(0, 0, 0, 0, 0, 4, 20, 26, 22, 10, 8, 3, 3, 5, 2, 0,0,0,0), len = 6:24) # 103 *
mar_1990 = data.frame('date' = as.POSIXct('1990-03-01'), 'val' = c(0, 0, 0, 0, 1, 1, 20, 14, 27, 14, 8, 4, 0, 0, 0, 0, 0, 0, 0), len = 6:24) # 89 *
apr_1990 = data.frame('date' = as.POSIXct('1990-04-01'), 'val' = c(0, 0, 0, 0, 0, 1, 6, 17, 17, 15, 14, 4, 4, 3, 0, 0, 0, 0, 0), len = 6:24) # 81 *
jun_1990 = data.frame('date' = as.POSIXct('1990-06-01'), 'val' = c(0, 0, 0, 0, 0, 0, 2, 13, 26, 19, 24, 13, 3, 3, 0, 1, 0, 0, 0), len = 6:24) # 104 *
aug_1990 = data.frame('date' = as.POSIXct('1990-08-01'), 'val' = c(0, 0, 0, 0, 0, 0, 1, 2, 6, 23, 26, 28, 9, 8, 2, 0, 0, 0, 0), len = 6:24) # 105 *
sep_1990 = data.frame('date' = as.POSIXct('1990-09-01'), 'val' = c(0, 0, 0, 0, 0, 0, 0, 1, 2, 5, 22, 27, 25, 7, 3, 4, 3, 1, 1), len = 6:24) # 101 *
oct_1990 = data.frame('date' = as.POSIXct('1990-10-01'), 'val' = c(0, 0, 0, 0, 1, 0, 0, 2, 6, 17, 17, 15, 7, 5, 5, 0, 0, 0, 0), len = 6:24) # 75
nov_1990 = data.frame('date' = as.POSIXct('1990-11-01'), 'val' = c(0, 0, 0, 0, 0, 1, 1, 2, 0, 3, 3, 8, 5, 2, 1, 0, 0, 0, 0), len = 6:24) # 26 * 
jan_1991 = data.frame('date' = as.POSIXct('1991-01-01'), 'val' = c(0, 0, 0, 0, 0, 0, 0, 3, 2, 0, 12, 32, 30, 8, 3, 0, 0, 0, 0), len = 6:24) # 90 * 
feb_1991 = data.frame('date' = as.POSIXct('1991-02-01'), 'val' = c(0, 0, 0, 0, 0, 0, 0, 0, 5, 7, 1, 6, 12, 14, 10, 0, 0, 0, 0), len = 6:24) # 55 *

## Combinging All of this data together into a single dataset
tagging_data = rbind(oct_1989, nov_1989, jan_1990, feb_1990, mar_1990, apr_1990, jun_1990, aug_1990, sep_1990, oct_1990, nov_1990, jan_1991, feb_1991)
## Making up pseudo tagging records. Columns are: date fish was caught, and length of the fish. 
pseudo_data = data.frame()
for(i in 1:length(tagging_data$date)){
  if(tagging_data$val[i] != 0){
    pseudo_data = rbind(pseudo_data, data.frame('date' = rep(tagging_data$date[i], times = tagging_data$val[i]), 'len' = rep(tagging_data$len[i], times = tagging_data$val[i])))
  }
}
colnames(pseudo_data) = c('date', 'length')

## Stripping out the month and year from each date object, then creating a new variable that is the month and year smushed together (matches histograms in moffit & parrish)
pseudo_data$month = months(pseudo_data$date)
pseudo_data$year = format(pseudo_data$date, "%Y")
pseudo_data$month_year = paste(pseudo_data$month, as.character(pseudo_data$year))

### Age at recruitment to juvenile fishing grounds. 
## Peak spawning occurs in July according to Luers, Demartini, Humphreys 2017
## Difference between first sampling trip (Oct) and peak spawning (July) = 3 months
age_at_recruitment = 3/12

## Estimating starting means from data
start_means = rbind(c(1, 17), c(11, 18), c(13, 20), c(14, 20), c(15, NA), c(15, NA), c(16, NA), c(17, NA), c(16, NA), c(10, 16), c(13, 16), c(13, 18), c(15, 17))
constrain_means = rbind(c(11, 17),c(11, 18), c(14, 19), c(13, 19), c(14, NA), c(14, NA), c(14, NA), c(17, NA), c(15, NA), c(10, 16.5), c(13, 17), c(13, 17), c(15, 19))

## Creating a function that we can use later for bootstrapping
length_freq_decomp = function(pseudo_data, plot = FALSE, fixed_modes = FALSE){
  lfdat = data.frame(stringsAsFactors = FALSE)
  for(i in 1:length(unique(pseudo_data$date))){
    print(i)
    curr_month_year = unique(pseudo_data$date)[i]
    month_data = pseudo_data[pseudo_data$date == curr_month_year, ]
    ## Determing the mean age of fish once they've recruited, assuming that they were born during peak spawning.
    mode.age = as.numeric(difftime(unique(month_data$date), as.POSIXct('1989-07-01'), "days")) / 365 # in years
    mode.age = c(mode.age, mode.age + 1) # Assumption is that if two cohorts are present, the second is one year older than the first
    # When we get to the second year of data, YOY for first year becomes 1+ year old, new cohort recruits. Because we based age on difftime for first cohort, we need to remove 1 year from all ages
    if(min(mode.age) > 1){ 
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
      if (format(month_data$date, "%m")[1] <= 7) { # if month is less than july (month of peak spawning) 
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
    lfdat$mode.se = as.numeric(levels(lfdat$mode.se)[lfdat$mode.se])
  }
  
  ## Sorting lfdat by the mean age. This makes it easier to visually inspect that fish are getting larger as they get older. This can get messed up during bimodal composition. 
  lfdat = lfdat[order(lfdat$mode.age), ]
  if(plot){
    plot(x = lfdat$mode.age, y = lfdat$mode.len, pch = 19)
  }
  return(lfdat)
}

### Creating table of fitted components for gausian and guassian mixture models
print('lfdat-ing')

lfdat = length_freq_decomp(pseudo_data, plot = TRUE, fixed_modes = TRUE)

##### Model Fitting #####
print('Fitting models')
results = list()
results$training = data.frame(stringsAsFactors = FALSE)
results$full = data.frame(stringsAsFactors = FALSE)

#### Fitting VB model
growth.ssnl.f<- growthvb.f
npf <- 1  #number of parameters passed to growth.ssnl.f (in this case k)
npA <- 2  #number of parameters in distribution of release ages for tag model

#### Fitting each data stream individually
# print('Estimating Model Parameters')
# ### 1. Mark Recapture Data
# 
# ## Specifying starting parameters, as well as upper and lower bounds for parameter estimation
# #        mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
# p0 <- c(  70,     1, .10,     1,   .10,      1,     0,   0,     1,     0)
# lb <- c(  50,   0.1, .05,   0.1,   .05,    0.1,     0,   0,     0,     0)
# ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,   0,     0,     0)
# 
# fit.vb.tagging.all.data <- nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=tagdat,otodat=otodat,lfdat=lfdat, wt.oto=0,wt.tag=1,wt.lf=0)
# results$full = rbind(results$full, cbind('Model 5 - Mark Recapture - All Data', t(as.vector(fit.vb.tagging.all.data$par))))
# 
# fit.vb.tagging.train.data <- nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=tagdat_train,otodat=otodat,lfdat=lfdat, wt.oto=0,wt.tag=1,wt.lf=0)$par
# results$training = rbind(results$training, cbind('Model 6 - Mark Recapture - Training Data', t(as.vector(fit.vb.tagging.train.data$par))))
# 
# ### 2. Length at Age Data
# 
# ## Setting intial params for otolith data
# #        mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
# p0 <- c(  70,     1, .10,     1,   .10,      1,     0,   0,     1,      1)
# lb <- c(  50,   0.1, .05,   0.1,   .05,    0.1,     0,  -10,  0.1,    0.1)
# ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,   10,   15,     15)
# 
# fit.vb.oto <- nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=tagdat_train,otodat=otodat,lfdat=lfdat, wt.oto=1,wt.tag=0,wt.lf=0)
# results$full = rbind(results$full, cbind('Length at Age', t(as.vector(fit.vb.oto$par))))
# 
# ### 3a. Length Frequency Data - Unconstrained Linf
# 
# ### First some notes about replicating Results of Moffitt and Parrish 1996 - ELEFAN model they used did not estimate a0. a0 is soaking up some of the observed variability that otherwise goes to K. In the function logl.lf.f within the script file joint_lkhd.r, uncommenting the line a0 = 0 will force this model.
# ### So which model is appropriate? Lets use AICc to find out
# 
# ## Unconstrained Fit
# ## Setting intial params for length frequency data
# #        mu.L, sig.L,  k, mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
# p0 <- c(  70,     1, .10,     1,   .10,      1,     0,   0,     1,      1)
# lb <- c(  50,   0.1, .05,   0.1,   .05,    0.1,     0,  -10,  0.1,    0.1)
# ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,   10,   15,     15)
# 
# fit.vb.lfu <- nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=tagdat_train,otodat=otodat,lfdat=lfdat, wt.oto=0, wt.tag=0, wt.lf=1)
# results$full = rbind(results$full, cbind('Length Frequency (Unconstrained)', t(as.vector(fit.vb.lfu$par))))
# 
# ### 3b. Length Frequency Data - Linf constrained by larger linf from oto/mr data
# 
# ## Constraining Linf to a constant - In this case, maximum Linf from oto or mark recapture
# if(is.factor(results$full[ ,2])){
#   lic = max(levels(results$full[ ,2])[results$full[ ,2]])
# } else {
#   lic = max(results$full[ ,2]) # Note: second column is mu.L parameter (mean of Linf)
# }
# # lic = 78 # Same as used by Moffitt and Parrish (1996)
# 
# ## Setting intial params for length frequency data
# #        mu.L, sig.L,  k, mu.A, sig.A, sig.sci, sig.f,   a0,  sig.oto, sig.lf
# p0 <- c(  lic,     1, .10,     1,   .10,      1,     0,   0,     1,      1)
# lb <- c(  lic,   0.1, .05,   0.1,   .05,    0.1,     0,  -10,  0.1,    0.1)
# ub <- c(  lic,  15.0, .50,   1.5,   .50,   15.0,     0,   10,   15,     15)
# 
# fit.vb.lfc <- nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=tagdat_train,otodat=otodat,lfdat=lfdat, wt.oto=0, wt.tag=0, wt.lf=1)
# results$full = rbind(results$full, cbind('Length Frequency (Constrained)', t(as.vector(fit.vb.lfc$par))))
# 
# ############## Is it appropriate to try to measure a0 using such limited data? #################
# ### Lets use AICc to find out
# ## AICc = 2k - 2log(L) + ((2k^2 + 2k) / (n-k-1))
# aicc_with_a0_and_sig.lf = 2*3 + 2*(40.02605) + ((2*3^2 + 2*3) / (21 - 3 - 1)) # 87.46386
# aicc_without_a0 = 2*2 + 2*(62.30009) + ((2*2^2 + 2*2) / (21 - 2 - 1)) # 129.2668
# aicc_without_sig.lf = 2*2 + 2*(359.5163) + ((2*2^2 + 2*2) / (21 - 2 - 1)) # 359.5163
# aicc_without_a0_or_sig.lf = 2*1 + 2*(7565.03) + ((2*1^2 + 2*1) / (21 - 1 - 1)) # 7565.03
# ### Conclusion: Yes, definitely, because AICc with a0 and sig.lf is more than 40 units lower
# ###############################################################################################
# 
# ### Setting intial params for all data
# #        mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
# p0 <- c(  70,     1, .10,     1,   .10,      1,     0,   0,     1,      1)
# lb <- c(  50,   0.1, .05,   0.1,   .05,    0.1,     0,  -10,  0.1,    0.1)
# ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,   10,   15,     15)
# 
# ### 7. Model including all Data sources - Equal weighting to each data type
# fit.vb.equalwt.grouped <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, otodat=otodat, lfdat=lfdat, wt.oto=1/dim(otodat)[1], wt.tag=1/dim(tagdat_train)[1], wt.lf=1/length(lfdat$curr_month_year))
# results$training = rbind(results$training, cbind('Model 7 - All Data - Equal Weighting', t(as.vector(fit.vb.equalwt.grouped$par))))
# fit.vb.equalwt.grouped <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, otodat=otodat, lfdat=lfdat, wt.oto=1/dim(otodat)[1], wt.tag=1/dim(tagdat)[1], wt.lf=1/length(lfdat$curr_month_year))
# results$full = rbind(results$full, cbind('Model 7 - All Data - Equal Weighting', t(as.vector(fit.vb.equalwt.grouped$par))))
# 
# ### 8. Model including all Data sources - weighting based on number of sample size
# fit.vb.byn.grouped <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, otodat=otodat, lfdat=lfdat, wt.oto=1, wt.tag=1, wt.lf=1)
# results$training = rbind(results$training, cbind('Model 8 - All Data - Weighted by n', t(as.vector(fit.vb.byn.grouped$par))))
# fit.vb.byn.grouped <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, otodat=otodat, lfdat=lfdat, wt.oto=1, wt.tag=1, wt.lf=1)
# results$full = rbind(results$full, cbind('Model 8 - All Data - Weighted by n', t(as.vector(fit.vb.byn.grouped$par))))
# 
# ### 9. Model including all Data sources treated individually - with equal weighting
# fit.vb.equalwt.indv <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 1/dim(otodat[otodat$source == 'ralston and miyamoto', ])[1], wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat_train)[1], wt.tag2 = 0, wt.lf = 1/length(lfdat$curr_month_year), wt.lf2 = 0)
# results$training = rbind(results$training, cbind('Model 9 - Separated Data - Equal Weighting', t(as.vector(fit.vb.equalwt.indv$par))))
# fit.vb.equalwt.indv <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 1/dim(otodat[otodat$source == 'ralston and miyamoto', ])[1], wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat)[1], wt.tag2 = 0, wt.lf = 1/length(lfdat$curr_month_year), wt.lf2 = 0)
# results$full = rbind(results$full, cbind('Model 9 - Separated Data - Equal Weighting', t(as.vector(fit.vb.equalwt.indv$par))))
# 
# ### 10. Model including all Data sources treated individually - weighting based on number of sample size
# fit.vb.byn.indv <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1, wt.oto2= 1, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)
# results$training = rbind(results$training, cbind('Model 10 - Separated Data - Weighted by n', t(as.vector(fit.vb.byn.indv$par))))
# fit.vb.byn.indv <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1, wt.oto2= 1, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)
# results$full = rbind(results$full, cbind('Model 10 - Separated Data - Weighted by n', t(as.vector(fit.vb.byn.indv$par))))
# 
# ### 11. Model without Ralston & Miyamoto - Equal weighting (Because Brett said this was shit!)
# fit.vb.byn.indv.no.ralston <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 0, wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat_train)[1], wt.tag2 = 0, wt.lf = 1/length(lfdat$curr_month_year), wt.lf2 = 0)
# results$training = rbind(results$training, cbind('Model 11 - Separated Data - Equal Weighting - No R&M', t(as.vector(fit.vb.byn.indv.no.ralston$par))))
# fit.vb.byn.indv.no.ralston <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 0, wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat)[1], wt.tag2 = 0, wt.lf = 1/length(lfdat$curr_month_year), wt.lf2 = 0)
# results$full = rbind(results$full, cbind('Model 11 - Separated Data - Equal Weighting - No R&M', t(as.vector(fit.vb.byn.indv.no.ralston$par))))
# 
# ### 12. Model without Ralston & Miyamoto - weighted by n (Because Brett said this was shit!)
# ub <- c( 80,  10, 3,   2,  .5,   15,     0,   5,   5,     5)
# 
# fit.vb.byn.indv.no.ralston <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1, wt.oto2= 0, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)
# results$training = rbind(results$training, cbind('Model 12 - Separated Data - Weighted by n - No R&M', t(as.vector(fit.vb.byn.indv.no.ralston$par))))
# fit.vb.byn.indv.no.ralston <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1, wt.oto2= 0, wt.oto3=1, wt.oto4 = 1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)
# results$full = rbind(results$full, cbind('Model 12 - Separated Data - Weighted by n - No R&M', t(as.vector(fit.vb.byn.indv.no.ralston$par))))
# 
# results$training$`time to 90%` = yrs_to_.9_linf(linf = as.numeric(levels(results$training$mu.L)[results$training$mu.L]), k = as.numeric(levels(results$training$k)[results$training$k]), a0 = as.numeric(levels(results$training$a0)[results$training$a0]))
# results$full$`time to 90%` = yrs_to_.9_linf(linf = as.numeric(levels(results$full$mu.L)[results$full$mu.L]), k = as.numeric(levels(results$full$k)[results$full$k]), a0 = as.numeric(levels(results$full$a0)[results$full$a0]))

#### Writing results out results to .csv  
# write.csv(results$training, file = file.path(run_results_dir, 'likelihood_parameter_estimates_with_training_data.csv'))
# write.csv(results$full, file = file.path(run_results_dir, 'likelihood_parameter_estimates_with_full_data.csv'))

#### Now we need to evaluate which model works the best
print('Evaluating Model Structures')
## We will do this by comparing each model's parameters from training data to observations in validation data
## Model scoring metric is as follows: sum((predicted - observed)^2) / n
## Lower scoring metric indicates better model fit

n_train    = round(dim(tagdat)[1] * (2/3))

evaluate_models = function(cross_validation_iterations = 10000){
mod_eval_results = data.frame(stringsAsFactors = FALSE)
model_na_results = rep(0, 7)
lit_vbgf = lit_vbgc_params[lit_vbgc_params$region %in% c('Hawaii - MHI & NWHI', 'Hawaii - MHI', 'Hawaii - NWHI'), ]
lit_vbgf_for_train = lit_vbgf[!(lit_vbgf$author %in% paste('Bayesian Model',1:4)), ]
bayes_models = lit_vbgf[(lit_vbgf$author %in% paste('Bayesian Model',1:4)), ]

mod_eval_results = foreach(i = 1:cross_validation_iterations, .combine = rbind) %dopar% {
  train_index = sample(1:dim(tagdat)[1], size = n_train, replace = FALSE)
  tagdat_train = tagdat[train_index, ]
  tagdat_validate = tagdat[-train_index, ]
  
  score = rep(NA, 7)
  ### Setting intial params for all data
  #        mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
  p0 <- c(  70,     1,  .10,     1,   .10,      1,     0,    0,      1,      1)
  lb <- c(  40,   0.01, .05,  0.01,   .05,   0.01,     0,  -15,   0.01,   0.01)
  ub <- c( 110,  15.0,  .50,   1.5,   .50,   15.0,     0,   15,     15,     15)
  
  var5 = NULL
  var5 = try(nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=tagdat_train, wt.oto=0,wt.tag=1,wt.lf=0)$par, silent = TRUE)
  if(is.numeric(var5[1])){
    score[1] = calculate_predictive_variance(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = var5[1], k = var5[3], Lr_obs = tagdat_validate[ ,2])
  } else {
    model_na_results[1] = model_na_results[1] + 1
  }
  
  var6 = NULL
    var6 = try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, otodat=otodat, lfdat=lfdat, wt.oto=1/dim(otodat)[1], wt.tag=1/dim(tagdat_train)[1], wt.lf=1/length(lfdat$curr_month_year))$par, silent = TRUE)
    if(is.numeric(var6[1])){
      score[2] = calculate_predictive_variance(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = var6[1], k = var6[3], Lr_obs = tagdat_validate[ ,2])
    } else {
      model_na_results[2] = model_na_results[2] + 1
    }

  var7 = NULL
    var7 = try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, otodat=otodat, lfdat=lfdat, wt.oto=1, wt.tag=1, wt.lf=1)$par, silent = TRUE)
    if(is.numeric(var7[1])){
      score[3] = calculate_predictive_variance(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = var7[1], k = var7[3], Lr_obs = tagdat_validate[ ,2])
    } else {
      model_na_results[3] = model_na_results[3] + 1
    }

  var8 = NULL
    var8 =  try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 1/dim(otodat[otodat$source == 'ralston and miyamoto', ])[1], wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat_train)[1], wt.tag2 = 0, wt.lf = 1/length(lfdat$curr_month_year), wt.lf2 = 0)$par, silent = TRUE)
      if(is.numeric(var8[1])){
      score[4] = calculate_predictive_variance(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = var8[1], k = var8[3], Lr_obs = tagdat_validate[ ,2])
    } else {
      model_na_results[4] = model_na_results[4] + 1
    }

  var9 = NULL
    var9 =  try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1, wt.oto2= 1, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)$par, silent = TRUE)
    if(is.numeric(var9[1])){
      score[5] = calculate_predictive_variance(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = var9[1], k = var9[3], Lr_obs = tagdat_validate[ ,2])
    } else {
      model_na_results[5] = model_na_results[5] + 1
    }

  var10 = NULL
    var10 =  try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 0, wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat_train)[1], wt.tag2 = 0, wt.lf = 1/length(lfdat$curr_month_year), wt.lf2 = 0)$par, silent = TRUE)
    if(is.numeric(var10[1])){
      score[6] = calculate_predictive_variance(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = var10[1], k = var10[3], Lr_obs = tagdat_validate[ ,2])
    } else {
      model_na_results[6] = model_na_results[6] + 1
    }

  var11 = NULL
    var11 = try(nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat = lfdat, lfdat2 = NULL, wt.oto = 1, wt.oto2= 0, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)$par, silent = TRUE)
    if(is.numeric(var11[1])){
      score[7] = calculate_predictive_variance(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = var11[1], k = var11[3], Lr_obs = tagdat_validate[ ,2])
    } else {
      model_na_results[7] = model_na_results[7] + 1
    }
    
  ## Now getting fits from literature data  
    ## lit_models
    lit_var_scores = rep(0, length(lit_vbgf_for_train$author))
    for(i in 1:length(lit_vbgf_for_train$author)){
      lit_var_scores[i] = calculate_predictive_variance(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = lit_vbgf_for_train$linf[i], k = lit_vbgf_for_train$k[i], Lr_obs = tagdat_validate[ ,2])
    }
   #  best_lit_score = min(lit_var_scores)
   #  best_lit_mod = lit_vbgf_for_train$author[which(lit_var_scores == min(lit_var_scores))]
    
    ## Bayes_models
    bayes_models = lit_vbgc_params[21:24, ]
    bayes_var_scores = rep(0, length(bayes_models$author))
    for(i in 1:length(bayes_models$author)){
      bayes_var_scores[i] = calculate_predictive_variance(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = bayes_models$linf[i], k = bayes_models$k[i], Lr_obs = tagdat_validate[ ,2])
    }
   # best_bayes_score = min(bayes_var_scores)
   #  best_bayes_mod = bayes_models$author[which(bayes_var_scores == min(bayes_var_scores))]
    
  ## Best Overall Model
  # best_model = c('Ensemble', 'Literature', 'Bayesian')[which.max(c(best_ll_score, best_lit_score, best_bayes_score))]
  new_line = c(score, lit_var_scores, bayes_var_scores)
  return(new_line)
  #return(data.frame(cbind(score, lit_var_scores, bayes_var_scores)))
  # return(data.frame('best_ll_mod' = best_ll_mod, 'best_lit_mod' = best_lit_mod, 'best_bayes_mod' = best_bayes_mod, 'll_score' = best_ll_score, 'lit_score' = best_lit_score, 'bayes_mod_score' = best_bayes_score, 'best_model' = best_model))
}
colnames(mod_eval_results) = c(paste('model', 5:11), lit_vbgf_for_train$author, bayes_models$author)
return(invisible(mod_eval_results))
}

#### What was prefered model?
n_iterations = 1000
mod_eval_results = as.data.frame(evaluate_models(cross_validation_iterations = n_iterations))
mod_eval_results_lf = as.data.frame(t(mod_eval_results[ ,1:7]))
mod_eval_results_lf$model_id = rownames(mod_eval_results_lf)
mod_eval_results_lf = reshape(mod_eval_results_lf, varying = colnames(mod_eval_results_lf[1:n_iterations]), idvar = 'model_id', direction = "long")
mod_eval_results_lf = mod_eval_results_lf

boxplot(mod_eval_results_lf$result ~ mod_eval_results_lf$model_id, ylim = c(0, 15))

#### Declaring the best model - The model that has the lowest mean evaluation result
## First finding best structure for integrative model
integrative_models = mod_eval_results[ ,2:7]
integrative_model_scores = c()
for(i in 1:dim(integrative_models)[1]){
  integrative_model_scores = c(integrative_model_scores, names(which.min(integrative_models[i, ])))
}
## Which model was most frequently the best one?
best_integrative_model = names(which.max(table(integrative_model_scores)))

## Now comparing best integrative model to just tagging data
integrative_vs_tagging = mod_eval_results[ ,which(colnames(mod_eval_results) %in% c('model 5', best_integrative_model))]
integrative_vs_tagging_model_scores = c()
for(i in 1:dim(integrative_models)[1]){
  integrative_vs_tagging_model_scores = c(integrative_vs_tagging_model_scores, names(which.min(integrative_vs_tagging[i, ])))
}
best_model = names(which.max(table(integrative_vs_tagging_model_scores)))

# 
# model_structure_selection = data.frame()
# nll_names = colnames(mod_eval_results)[1:7]
# lit_names = colnames(mod_eval_results)[8:18]
# bayes_names = colnames(mod_eval_results)[19:22]
# for(i in 1:length(mod_eval_results[ ,1])){
#   score_ens = min(mod_eval_results[i,1:7], na.rm = TRUE)
#   best_ens = nll_names[which(mod_eval_results[i,1:7] == score_ens)]
#   score_lit = min(mod_eval_results[i,8:18], na.rm = TRUE)
#   best_lit = lit_names[which(mod_eval_results[i,8:18] == score_lit)]
#   score_bayes = min(mod_eval_results[i,19:22], na.rm = TRUE)
#   best_bayes = bayes_names[which(mod_eval_results[i,19:22] == score_bayes)]
#   best_overall = c('MLE', 'Lit', 'Bayes')[which.min(c(score_ens, score_lit, score_bayes))]
#   best_mod = c(best_ens, best_lit, best_bayes)[which.min(c(score_ens, score_lit, score_bayes))]
#   write_line = data.frame('best_ll_mod' = best_ens, 'score_ensemble' = score_ens, 'best_lit_mod' = best_lit, 'score_lit' = score_lit, 'best_bayes_mod' = best_bayes, 'score_bayes' = score_bayes, 'best_model' = best_overall, 'best_indv_mod' = best_mod)
#   model_structure_selection = rbind(model_structure_selection, write_line)
# }
# 
# mod_eval_results_table = aggregate(model_structure_selection$best_ll_mod, by = list(model_structure_selection$best_ll_mod), FUN = length)
# colnames(mod_eval_results_table) = c('model', 'n')
# 
# all_eval_results_table = aggregate(model_structure_selection$best_indv_mod, by = list(model_structure_selection$best_indv_mod), FUN = length)
# colnames(all_eval_results_table) = c('model', 'n')

send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste('best MLE model:', mod_eval_results_table$model[which.max(mod_eval_results_table$n)]))
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste('best model overall:', all_eval_results_table$model[which.max(all_eval_results_table$n)]))
# 
# ens.results = as.vector(mod_eval_results[ ,1:7])[!is.na(as.vector(mod_eval_results[ ,1:7]))]
# ens.mod.s.range = range(ens.results)
# ens.mod.s.mean = mean(ens.results)
# ens.mod.s.se = std_error(ens.results)
# 
# lit.results = as.vector(mod_eval_results[ ,8:18])[!is.na(as.vector(mod_eval_results[ ,8:18]))]
# lit.mod.s.range = range(lit.results)
# lit.mod.s.mean = mean(lit.results)
# lit.mod.s.se = std_error(lit.results)
# 
# bayes.results = as.vector(mod_eval_results[ ,19:22])[!is.na(as.vector(mod_eval_results[ ,19:22]))]
# bayes.mod.s.range = range(bayes.results)
# bayes.mod.s.mean = mean(bayes.results)
# bayes.mod.s.se = std_error(bayes.results)

# ## Plotting this
# pdf(file.path(run_results_dir, 'Barplot of Best Model Structures fit to Test Data.pdf'), width = 11, height = 8.5)
#   barplot(prop.table(table(model_structure_selection$best_ll_mod)))
# dev.off()
# 
# pdf(file.path(run_results_dir, 'Barplot Comparing Bayes model fits to Test Data.pdf'), width = 11, height = 8.5)
#   barplot(prop.table(table(model_structure_selection$best_bayes_mod)))
# dev.off()
# 
#   pdf(file.path(run_results_dir, 'Barplot Comparing Lit fits to Test Data.pdf'), width = 11, height = 8.5)
#   barplot(prop.table(table(model_structure_selection$best_lit_mod)))
# dev.off()
# 
#   pdf(file.path(run_results_dir, 'Barplot Comparing NLL Bayes and Lit fits to Test Data.pdf'), width = 11, height = 8.5)
#   barplot(prop.table(table(model_structure_selection$best_model)))
# dev.off()

## Write results out
save.image(file = file.path(run_results_dir, 'workspace_image_preboot.RData'))

### Now bootstrapping the best model
print('Bootstrapping prefered model')
#send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = 'Booting Prefered Model')
boot_iterations = 10000
bootstrap_results = list()

mod_eval_results_table = aggregate(model_structure_selection$best_ll_mod, by = list(model_structure_selection$best_ll_mod), FUN = length)

#### Bootstrap out tagging only model
## Setting intial params for all data
#        mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
p0 <- c(  70,     1, .10,     1,    .1,      1,     0,   0,    1,      1)
lb <- c(  50,   0.1, .05,   0.1,   .05,    0.1,     0, -10,  0.1,    0.1)
ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,  10,   15,     15)

mod_eval_results_table[which.max(mod_eval_results_table[ ,2]), 1] == 'model 5'
  ## Specifying starting parameters, as well as upper and lower bounds for parameter estimation
  #        mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
  p0 <- c(  70,     1, .10,   1.0,   .10,      1,     0,   0,    0,     0)
  lb <- c(  50,   0.1, .05,   0.1,   .05,    0.1,     0,   0,    0,     0)
  ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,   0,    0,     0)
  
  print('Booting Model 5')
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 5")
  timer6full = proc.time()
  bootstrap_results$booted_param_ests_all_tagging_data = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_5', boot_iterations = boot_iterations, wt.oto = 0, wt.lf = 0, wt.tag = 1, tagdat = tagdat)
  boot_time =  (proc.time() -  timer6full)[3] / 60 / 60
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 5  complete!"))
  
### Now bootstrapping our prefered model  
if (best_model == 'model 6') {
  ## 6. Model including all Data sources - Equal weighting to each data type
  print('Booting Model 6')
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 6")
  timer7 = proc.time()
  bootstrap_results$booted_param_ests_model6 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_6',boot_iterations = boot_iterations, wt.oto = 1/length(otodat$age), wt.lf = 1/length(lfdat$curr_month_year), wt.tag = 1/dim(tagdat)[1], otodat = otodat, tagdat = tagdat, pseudolf = pseudo_data)
  boot_time =  (proc.time() -  timer7)[3] / 60 / 60
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 6 complete!"))
  
} else if (best_model == 'model 7') {
  ## 7. Model including all Data sources - weighting based on number of sample size
  print('Booting Model 7')
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 7")
  timer8 = proc.time()
  bootstrap_results$booted_param_ests_model7 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_7_all_data', boot_iterations = boot_iterations,tagdat=tagdat, otodat=otodat, pseudolf=pseudo_data, wt.oto=1, wt.tag=1, wt.lf=1)
  boot_time =  (proc.time() -  timer8)[3] / 60 / 60
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 7 complete!"))
  
} else if (best_model == 'model 8') {
  ## 8. Model including all Data sources treated individually - with equal weighting
  print('Booting Model 8')
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 8")
  timer9 = proc.time()
  bootstrap_results$booted_param_ests_model8 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_8_all_data', boot_iterations = boot_iterations, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], pseudolf=pseudo_data, pseudolf2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 1/dim(otodat[otodat$source == 'ralston and miyamoto', ])[1], wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat)[1], wt.tag2 = 0, wt.lf = 1/length(pseudolf$curr_month_year), wt.lf2 = 0)
  boot_time =  (proc.time() -  timer9)[3] / 60 / 60
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 8 complete!"))
  
} else if (best_model == 'model 9') {
  ## 9. Model including all Data sources treated individually - weighting based on number of sample size
  print('Booting Model 9')
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 9")
  timer10 = proc.time()
  bootstrap_results$booted_param_ests_model9 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_9', boot_iterations = boot_iterations, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], pseudolf=pseudo_data, pseudolf2=NULL, wt.oto= 1, wt.oto2= 1, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)
  boot_time =  (proc.time() -  timer10)[3] / 60 / 60
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 9 complete!"))
  
} else if (best_model == 'model 10') {
  ## 10. Model without Ralston & Miyamoto - Equal weighting (Because Brett said this was shit!)
  print('Booting Model 10')
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 10")
  timer11 = proc.time()
  bootstrap_results$booted_param_ests_model10 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_10', boot_iterations = boot_iterations, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], pseudolf=pseudo_data, pseudolf2 = NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 0, wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat)[1], wt.tag2 = 0, wt.lf = 1/length(pseudolf$curr_month_year), wt.lf2 = 0)
  boot_time =  (proc.time() -  timer11)[3] / 60 / 60
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 10 complete!"))
  
} else if (best_model == 'model 11') {
  ### 11. Model without Ralston & Miyamoto - weighted by n (Because Brett said this was shit!)
  print('Booting Model 11')
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 11")
  timer12 = proc.time()
  bootstrap_results$booted_param_ests_model11 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_11', boot_iterations = boot_iterations, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], pseudolf=pseudo_data, pseudolf2=NULL, wt.oto= 1, wt.oto2= 0, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)
  boot_time =  (proc.time() -  timer12)[3] / 60 / 60
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 11 complete!"))
}

#save.image(file = file.path(run_results_dir, 'workspace_image.RData'))
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = 'Run complete!')

#send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "If you get this before boot complete, something is fucked!"))


# 
# mod_struc_selection_pred = data.frame(NULL)
# mod_struc_selection_resid = data.frame(NULL)
# for(i in 1:length(results$training$Data)){
#   l_inf = as.numeric(levels(results$training$mu.L)[i])
#   k =  as.numeric(levels(results$training$k)[i])
#   a0 = as.numeric(levels(results$training$a0)[i])
#   mod_struc_selection_pred = rbind(mod_struc_selection_pred, predict_recapture_length(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], l_inf, k = k, a0))
#   mod_struc_selection_resid = rbind(mod_struc_selection_resid, mod_struc_selection_pred[i, ] - tagdat_validate[ ,2])
# }
# rownames(mod_struc_selection_pred)  = results$training$Data
# rownames(mod_struc_selection_resid) = results$training$Data
# training_mod_scores = cbind(results$training$Data, rowSums(mod_struc_selection_resid^2)/dim(mod_struc_selection_resid)[2])
# 
# #### Lets turn this into a pretty plot
# pdf(file.path(run_results_dir, 'Predicted vs. Observed LR with validation data for models 6-12.pdf'), width = 11, height = 8.5)
# par(mfrow = c(2, ceiling(length(rownames(mod_struc_selection_pred))/2)))
# for(i in 1:nrow(mod_struc_selection_pred)){
#   model_id = strsplit(rownames(mod_struc_selection_pred)[i], split = '-')[[1]][1]
#   plot(y = mod_struc_selection_pred[i, ], x = tagdat_validate[ ,2],
#        xlab = 'Observed Length at Recapture (cm)', xlim = c(15, 80), 
#        ylab = 'Predicted Length at Recapture (cm)', ylim = c(15, 80),
#        main = paste(model_id,'\n Linf = ', round(as.numeric(levels(results$training$mu.L)[i]), digits = 2),' k = ', round(as.numeric(levels(results$training$k)[i]), digits = 2), sep = ""))
#   abline(1, 1, lty = 2)
#   arrows(x0 = 60, y0 = 70, x1 = 65, y1 = 70, length = 0.1, angle = 30, code = 2)
#   text(x = 45, y = 70, labels = "Line of 1:1 agreement", cex = .75)
#   model_var = sum(mod_struc_selection_resid[i, ]^2)/ length(tagdat_validate[ ,2])
#   text(x = 60, y = 30, labels = paste("Predictive Variance:", round(model_var, digits = 3)), cex = .75)
# }
# dev.off()

# 
# #### Rerunning best fit models (8 and 10) with full data
# ### 8. Model including all Data sources - weighting based on number of sample size
# fit.vb.byn.grouped_full <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, otodat=otodat, lfdat=lfdat, wt.oto=1, wt.tag=1, wt.lf=1)
# results$full = rbind(results$full, cbind('Model 8 - All Data - Weighted by n', t(as.vector(fit.vb.byn.grouped$par))))
# 
# ### 10. Model including all Data sources treated individually - weighting based on number of sample size
# fit.vb.byn.indv <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1, wt.oto2= 1, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)
# results$full = rbind(results$full, cbind('Model 10 - Separated Data - Weighted by n', t(as.vector(fit.vb.byn.indv$par))))
# 
# colnames(results$training) = c('Data',  'mu.L', 'sig.L',  'k',  'mu.A', 'sig.A', 'sig.sci', 'sig.f', 'a0', 'sig.oto', 'sig.lf' )
# colnames(results$full) = c('Data',  'mu.L', 'sig.L',  'k',  'mu.A', 'sig.A', 'sig.sci', 'sig.f', 'a0', 'sig.oto', 'sig.lf' )

# 
# 
# ### Now we specify the final number of successfully bootstrapped estimates we want to obtained
# boot_iterations = 10000
# 
# ### Now we bootstrap our tagging data only model (Model 1)
# bootstrap_results = list()
# 
# ## Specifying starting parameters, as well as upper and lower bounds for parameter estimation
# #        mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
# p0 <- c(  70,     1, .10,   1.0,   .10,      1,     0,   0,    0,     0)
# lb <- c(  50,   0.1, .05,   0.1,   .05,    0.1,     0,   0,    0,     0)
# ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,   0,    0,     0)
# 
# ## Bootstrapping parameter estimates from full tagging data set
# 
# ## Bootstrapping parameter estimates from training tagging data set
# print('Booting Model 1 - Training Data')
# timer1atrain = proc.time()
# bootstrap_results$booted_param_ests_model1_train = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_1_training_data', boot_iterations = boot_iterations, wt.oto = 0, wt.lf = 0, wt.tag = 1, tagdat = tagdat_train)
# boot_time =  (proc.time() -  timer1atrain)[3] / 60 / 60
# #send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 1 - Training Data complete!"))
# 
# ### Now we bootstrap our combined models (Model 7, b, c, d)
# ## Setting intial params for all data
# #        mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
# p0 <- c(  70,     1, .10,     1,    .1,      1,     0,   0,    1,      1)
# lb <- c(  50,   0.1, .05,   0.1,   .05,    0.1,     0, -10,  0.1,    0.1)
# ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,  10,   15,     15)
# 
# ### Now we bootstrap our combined models on training data
# 
#  
#    
#   
#   
#   ##### So after we did this, we compared notes on models and found the ones that were the most bestest.
#           ## Those were Models 8 and 10.
#   #### So we bootstrap them again but on the full data set
#   results_all_data = data.frame()
#   results_all_data = rbind(results_all_data, cbind('Mark Recapture - Training Data', t(as.vector(fit.vb.tagging.all.data$par))))
#   
#   ## Setting intial params for all data
#   #        mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
#   p0 <- c(  70,     1, .10,     1,    .1,      1,     0,   0,    1,      1)
#   lb <- c(  50,   0.1, .05,   0.1,   .05,    0.1,     0, -10,  0.1,    0.1)
#   ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,  10,   15,     15)
#   
#   ### 8. Model including all Data sources - weighting based on number of sample size
#   fit.vb.byn.grouped <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, otodat=otodat, pseudo_data=pseudo_data, wt.oto=1, wt.tag=1, wt.lf=1)
#   results_all_data = rbind(results_all_data, cbind('All Data - Weighted by n', t(as.vector(fit.vb.byn.grouped$par))))
#   
#   print('Booting Model 8')
#   timer8 = proc.time()
#   bootstrap_results$booted_param_ests_model8_all_data = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_8_all_data', boot_iterations = boot_iterations, wt.oto = 1, wt.lf = 1, wt.tag = 1, otodat = otodat, tagdat = tagdat, pseudolf = pseudo_data)
#   boot_time =  (proc.time() -  timer8)[3] / 60 / 60
#   #send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 8 complete!"))
#   
#   
#   ### 9. Model including all Data sources treated individually - with equal weighting
#   fit.vb.equalwt.indv <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat_train, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 1/dim(otodat[otodat$source == 'ralston and miyamoto', ])[1], wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat_train)[1], wt.tag2 = 0, wt.lf = 1/length(lfdat$curr_month_year), wt.lf2 = 0)
#   results_all_data = rbind(results_all_data, cbind('Separated Data - Equal Weighting', t(as.vector(fit.vb.equalwt.indv$par))))
#   
#   print('Booting Model 10')
#   timer9 = proc.time()
#   bootstrap_results$booted_param_ests_model9_all_data = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_10_all_data', boot_iterations = boot_iterations, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], pseudolf=pseudo_data, lfdat2=NULL, wt.oto= 1, wt.oto2= 1, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)
#   boot_time =  (proc.time() -  timer9)[3] / 60 / 60
#   #send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 10 complete!"))
#   
#   
#   
#   
  
 #### Once everything has been booted, we need to compare parameter estimates

## Plotting predicted vs. observed length at recapture
### Constructing data frame of predicted values and residuals under each model
predicted_recapture_lengths = data.frame(NULL)
predicted_recapture_residuals = data.frame()
for(i in 1:length(names(bootstrap_results))){
  med_boot_l = bootstrap_results[[i]]$boot_stats['Median','mu.L']
  med_boot_k = bootstrap_results[[i]]$boot_stats['Median','k']
  med_boot_a0 = bootstrap_results[[i]]$boot_stats['Median','mu.A']
  predicted_recapture_lengths = rbind(predicted_recapture_lengths, predict_recapture_length(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = bootstrap_results[[i]]$boot_stats['Median','mu.L'], k = bootstrap_results[[i]]$boot_stats['Median','k'], a = bootstrap_results[[i]]$boot_stats['Median','mu.A']))
  predicted_recapture_residuals = rbind(predicted_recapture_residuals, predict_recapture_length(Lm = tagdat_validate[ ,1], dt = tagdat_validate[ ,4], linf = bootstrap_results[[i]]$boot_stats['Median','mu.L'], k = bootstrap_results[[i]]$boot_stats['Median','k']) - tagdat_validate[ ,2])
}
rownames(predicted_recapture_lengths)   = names(bootstrap_results)
rownames(predicted_recapture_residuals) = names(bootstrap_results)

## Plotting results
pdf(file.path(run_results_dir, 'Predicted vs. Observed LR with validation data.pdf'), width = 8.5, height = 11)
  par(mfrow = c(length(rownames(predicted_recapture_lengths))/2, 2))
  for(i in 1:nrow(predicted_recapture_lengths)){
    name_decomp = strsplit(names(bootstrap_results)[i], split = 'model')[[1]]
    model_name = paste('model', name_decomp[length(name_decomp)])
    plot(y = predicted_recapture_lengths[i, ], x= tagdat_validate[ ,2],
         xlab = 'Observed Length at Recapture (cm)', xlim = c(15, 80),
         ylab = 'Predicted Length at Recapture (cm)', ylim = c(15, 80),
         main = paste(model_name, '\n Linf = ', round(bootstrap_results[[i]]$boot_stats['Median','mu.L'], digits = 2),'k = ', round(bootstrap_results[[i]]$boot_stats['Median','k'], digits = 2), sep = ""))
    abline(1, 1, lty = 2)
    arrows(x0 = 60, y0 = 70, x1 = 65, y1 = 70, length = 0.1, angle = 30, code = 2)
    text(x = 45, y = 70, labels = "Line of 1:1 agreement", cex = .75)
    text(x = 50, y = 30, labels = paste("Variance Coefficient:", round(var(x = as.numeric(predicted_recapture_lengths[i, ]), y = as.numeric(tagdat_validate[ ,2])), digits = 3)), cex = .75)
  }
dev.off()


### Plotting predicted vs. observed residuals
pdf(file.path(run_results_dir, 'Residual vs. Observed LR with validation data.pdf'), width = 8.5, height = 11)
  par(mfrow = c(length(rownames(predicted_recapture_lengths))/2, 2))
  for(i in 1:nrow(predicted_recapture_lengths)){
    name_decomp = strsplit(names(bootstrap_results)[i], split = 'model')[[1]]
    model_name = paste('model', name_decomp[length(name_decomp)])
    plot(as.numeric(predicted_recapture_residuals[i, ]) ~ tagdat_validate[ ,2],
         xlab = 'Observed Length at Recapture (cm)', xlim = c(15, 80),
         ylab = 'Residual Length at Recapture (cm)', ylim = c(-1*max(abs(range(predicted_recapture_residuals))), max(abs(range(predicted_recapture_residuals)))),
         main = paste(model_name, '\n Linf = ', round(bootstrap_results[[i]]$boot_stats['Median','mu.L'], digits = 2),'k = ', round(bootstrap_results[[i]]$boot_stats['Median','k'], digits = 2), sep = ""))
    abline(0, 0, lty = 2)
  }
dev.off()


### Plotting each parameter estimate
for(r in 1:length(colnames(bootstrap_results$booted_param_ests_model7$boot_stats))){
  pdf(file.path(run_results_dir, paste('Parameter Estimate Boxplots-',colnames(bootstrap_results$booted_param_ests_model7$boot_stats)[r],'.pdf', sep = "")), height = 8.5, width = 11)
  par(mfrow = c(2, 1))
  lboxdata = list()
  for(i in 1:length(names(bootstrap_results))){
    lboxdata[[i]] = bootstrap_results[[i]]$raw_boot_data[, r][!is.na(bootstrap_results[[i]]$raw_boot_data[, r])]
  }
  names(lboxdata) = names(bootstrap_results)
  boxplot(lboxdata, main = colnames(bootstrap_results$booted_param_ests_model7$boot_stats)[r])
  dev.off()
  }

#### Script wrapup
save.image(file = file.path(run_results_dir, 'workspace_image.RData'))
script_time = round((proc.time() - script_timer)[3]/60/60, digits = 3)
print(paste('Total Time:', script_time))
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(script_time, "Hours later, All Done!"))






# #### Confirmation of Length Frequency Results since they're not what's published
# library('TropFishR')
# data("synLFQ7")
# x = synLFQ7
# x$sample.no = 1:13
# x$midLengths = seq(8,24)
# x$dates = unique(as.Date(pseudo_data$date))
# 
# catch = aggregate(pseudo_data$year, by = list(pseudo_data$date, pseudo_data$length), FUN = length)
# x$catch = as.matrix(xtabs(x~Group.2+Group.1, data=catch))
# colnames(x$catch) = NULL
# rownames(x$catch) = NULL
# 
# synLFQ7a <- lfqModify(x, bin_size = 4)
# lfqbin <- lfqRestructure(synLFQ7a, MA = 3, addl.sqrt = FALSE)
# opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
# plot(lfqbin, Fname = "catch", date.axis = "modern")
# plot(lfqbin, Fname = "rcounts", date.axis = "modern")
# par(opar)
# 
# res_KScan <- ELEFAN(synLFQ7a, Linf_fix = 68,
#                     MA=5, addl.sqrt = FALSE, hide.progressbar = TRUE)
# 
# res_RSA <- ELEFAN(synLFQ7a, Linf_range = seq(60,70,1), MA = 5,
#                   K_range = seq(0.01,2,0.1), addl.sqrt = TRUE,
#                   hide.progressbar = TRUE, contour=5)
# 
# n <- length(res_RSA$score_mat)
# best_scores <- sort(res_RSA$score_mat,partial=n-0:2)[n-0:2]
# ind <- arrayInd(which(res_RSA$score_mat %in% best_scores),
#                 dim(res_RSA$score_mat))
# Ks <- as.numeric(rownames(res_RSA$score_mat)[ind[,1]])
# Linfs <- as.numeric(colnames(res_RSA$score_mat)[ind[,2]])
# 
# res_loop <- vector("list", 3)
# for(i in 1:3){
#   tmp <- ELEFAN(synLFQ7a,
#                 Linf_range = seq(Linfs[i]-2, Linfs[i]+2, 0.2),
#                 K_range = seq(Ks[i]-0.1, Ks[i]+0.1, 0.05),
#                 MA = 5,
#                 addl.sqrt = TRUE,''';;;klHk;
#                 hide.progressbar = TRUE,
#                 contour=5)
#   res_loop[[i]] <- cbind(Rn_max=tmp$Rn_max, t(as.matrix(tmp$par)))
# }
# results <- do.call(rbind, res_loop)
# 
# r = ELEFAN(x, Linf_fix = 100)

## Debug parameters
#otodat = otodat[otodat$source == 'demartini', ]; wt.oto = 1; otodat2 = otodat[otodat$source == 'ralston and miyamoto', ]; wt.oto2 = 0; otodat3 = otodat[otodat$source == 'andrews bomb carbon', ]; wt.oto3 = 1; otodat4 = otodat[otodat$source == 'andrews lead radium', ]; wt.oto4 = 1; pseudolf = pseudo_data; wt.lf = 1; tagdat = tagdat_train; wt.tag = 1; tagdat2 = NULL; wt.tag2 = NULL; pseudolf2 = NULL; wt.lf2 = NULL

# 
# schnute_likelihood = function(a, b, c, alpha){
#   linf * (1 + alpha * exp(-a*t^c))
# }
