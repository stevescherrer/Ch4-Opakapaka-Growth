# HPCC

## Declaring Directory Path
proj_dir = getwd()
print(proj_dir)
# if(!"Okamoto_Mark_Recapture" %in% strsplit(proj_dir, '/')[[1]]){
#   proj_dir = file.path(getwd(), "Okamoto_Mark_Recapture")
# }
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
#### Mark Recapture Data 
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

n_recaps = data.frame('recapture events' = unique(paka_growth$n_recaptures), 'n_fish' = 0)
i = 1
for(n_recap in unique(paka_growth$n_recaptures)){
  n_recaps$n_fish[i] = dim(paka_growth[paka_growth$n_recaptures == n_recap, ])[1]
  i = i+1
}

#### Creating a subset data frame that removes recording errors in length and time
# paka_growth = subset(paka_growth, dL > 0)
length(which(paka_growth$dt < 60/365)) #46
paka_growth = subset(paka_growth, dt >= 60/365)
tagdat = as.matrix(data.frame('L1' = paka_growth$Lm, "L2" = paka_growth$Lr, " " = rep(0, length(paka_growth$Lr)), "dt" = paka_growth$dt, "L2measurer" = rep(0, length(paka_growth$Lr))))

#### Creating Second tagging dataset from PIFG data
pifg20072013 = read.csv(file.path(data_dir, 'PIFG 2007-2013.csv'), stringsAsFactors = FALSE)
pifg20072013$rel_date = as.POSIXct(pifg20072013$rel_date, format = "%m/%d/%Y")
pifg20072013$recap_date = as.POSIXct(pifg20072013$recap_date, format = "%m/%d/%Y")
pifg20072013$dt = difftime(pifg20072013$recap_date, pifg20072013$rel_date)

### 2014-2015 data
pifg20142015 = read.csv(file.path(data_dir, 'PIFG 2014-2015.csv'), stringsAsFactors = FALSE)
pifg20142015$rel_date = as.POSIXct(pifg20142015$rel_date, format = "%m/%d/%Y")
pifg20142015$recap_date = as.POSIXct(pifg20142015$recap_date, format = "%m/%d/%Y")
pifg20142015$rel_FL[pifg20142015$Units == 'in'] = pifg20142015$rel_FL[pifg20142015$Units == 'in'] * in_to_cm
pifg20142015$recap_FL[pifg20142015$Units == 'in'] = pifg20142015$recap_FL[pifg20142015$Units == 'in'] * in_to_cm
pifg20142015$dt = difftime(pifg20142015$recap_date, pifg20142015$rel_date)

pifg_data = data.frame('L1' = c(pifg20072013$rel_FL, pifg20142015$rel_FL), 'L2' = c(pifg20072013$recap_FL, pifg20142015$recap_FL), " " = rep(0, length(c(pifg20072013$recap_FL, pifg20142015$recap_FL))), 'dt' = c(pifg20072013$dt, pifg20142015$dt) / 365, "L2measurer" = rep(0, length(c(pifg20072013$recap_FL, pifg20142015$recap_FL))))

## Removing any pifg data with time at liberty < 60 days
pifg_data = pifg_data[pifg_data$dt >= 60/365, ]

tagdat2 = pifg_data


#### Otolith Data (Ralston and Miyamoto 1983, DeMartini et al. 1994, Andrews et al. 2012) 
otodat = read.csv(file.path(data_dir, "RalstonMiyamotoandDemartiniAndrews.csv"))
colnames(otodat) = c("age", "len", "source")

#### Length Frequency Data 
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
results = data.frame(stringsAsFactors = FALSE)

#### Fitting Full models without bootstrapping
growth.ssnl.f<- growthvb.f
npf <- 1  #number of parameters passed to growth.ssnl.f (in this case k)
npA <- 2  #number of parameters in distribution of release ages for tag model


#### Fitting each data stream individually
print('Estimating Model Parameters')
### 1. Mark Recapture Data
#Specifying starting parameters, as well as upper and lower bounds for parameter estimation
##       mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
p0 <- c(  70,     1, .10,     1,   .10,      1,     0,   0,     1,     0)
lb <- c(  50,   0.1, .05,   0.1,   .05,    0.1,     0,   0,     0,     0)
ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,   0,     0,     0)

fit.vb.tagging.all.data <- nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=tagdat,otodat=otodat,lfdat=lfdat, wt.oto=0,wt.tag=1,wt.lf=0)
results = rbind(results, cbind('Model 5 - Mark Recapture - All Data', t(as.vector(fit.vb.tagging.all.data$par))))

### Setting intial params for all data
##       mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
p0 <- c(  70,     1, .10,     1,   .10,      1,     0,   0,     1,      1)
lb <- c(  40,   0.01, .05,   0.1,   .05,    0.1,     0,  -10,  0.1,    0.1)
ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,   10,   15,     15)

### 6. Model including all Data sources - Equal weighting to each data type
fit.vb.equalwt.grouped <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, otodat=otodat, lfdat=lfdat, wt.oto=1/dim(otodat)[1], wt.tag=1/dim(tagdat)[1], wt.lf=1/length(lfdat$curr_month_year))
results = rbind(results, cbind('Model 6 - All Data - Equal Weighting', t(as.vector(fit.vb.equalwt.grouped$par))))

### 7. Model including all Data sources - weighting based on number of sample size
fit.vb.byn.grouped <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, otodat=otodat, lfdat=lfdat, wt.oto=1, wt.tag=1, wt.lf=1)
results = rbind(results, cbind('Model 7 - All Data - Weighted by n', t(as.vector(fit.vb.byn.grouped$par))))

### 8. Model including all Data sources treated individually - with equal weighting
fit.vb.equalwt.indv <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 1/dim(otodat[otodat$source == 'ralston and miyamoto', ])[1], wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat)[1], wt.tag2 = 0, wt.lf = 1/length(lfdat$curr_month_year), wt.lf2 = 0)
results = rbind(results, cbind('Model 8 - Separated Data - Equal Weighting', t(as.vector(fit.vb.equalwt.indv$par))))

### 9. Model including all Data sources treated individually - weighting based on number of sample size
fit.vb.byn.indv <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1, wt.oto2= 1, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)
results = rbind(results, cbind('Model 9 - Separated Data - Weighted by n', t(as.vector(fit.vb.byn.indv$par))))

### 10. Model without Ralston & Miyamoto - Equal weighting (Because Brett said this was shit!)
fit.vb.byn.indv.no.ralston <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 0, wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat)[1], wt.tag2 = 0, wt.lf = 1/length(lfdat$curr_month_year), wt.lf2 = 0)
results = rbind(results, cbind('Model 10 - Separated Data - Equal Weighting - No R&M', t(as.vector(fit.vb.byn.indv.no.ralston$par))))

### 11. Model without Ralston & Miyamoto - weighted by n (Because Brett said this was shit!)
fit.vb.byn.indv.no.ralston <- nlminb(p0, joint.logl.f, lower=lb, upper=ub, npf=npf, npA=npA, tagdat=tagdat, tagdat2 = NULL, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], lfdat=lfdat, lfdat2=NULL, wt.oto= 1, wt.oto2= 0, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 0, wt.lf = 1, wt.lf2 = 0)
results = rbind(results, cbind('Model 11 - Separated Data - Weighted by n - No R&M', t(as.vector(fit.vb.byn.indv.no.ralston$par))))


##### Bootstrapping alternative models #### (not preffered)
bootstrap_results = list()
n_iterations = 10000

### Now we'll bootstrap the prefered model structure  
## Setting intial params for all data
#        mu.L, sig.L,  k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
p0 <- c(  70,     1, .10,     1,    .1,      1,     0,   0,    1,      1)
lb <- c(  50,   0.1, .05,   0.1,   .05,    0.1,     0, -10,  0.1,    0.1)
ub <- c( 110,  15.0, .50,   1.5,   .50,   15.0,     0,  10,   15,     15)

print('Booting Model 6')
timer6 = proc.time()
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 6")
bootstrap_results$booted_param_ests_model6 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_6',boot_iterations = boot_iterations, wt.oto = 1/length(otodat$age), wt.lf = 1/length(lfdat$curr_month_year), wt.tag = 1/dim(tagdat)[1],   otodat = otodat, tagdat = tagdat, pseudolf = pseudo_data)
boot_time =  (proc.time() -  timer6)[3] / 60 / 60
send_push(user = 'uGEHvA4hr36tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 6 complete!"))

print('Booting Model 7')
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 7")
timer7 = proc.time()
bootstrap_results$booted_param_ests_model7 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_7_all_data', boot_iterations = boot_iterations,tagdat=tagdat, otodat=otodat, pseudolf=pseudo_data, wt.oto=1, wt.tag=1, wt.lf=1)
boot_time =  (proc.time() -  timer7)[3] / 60 / 60
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 7 complete!"))

print('Booting Model 8')
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 8")
timer8 = proc.time()
bootstrap_results$booted_param_ests_model8 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_8_all_data', boot_iterations = boot_iterations, tagdat=tagdat,  otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], pseudolf=pseudo_data, pseudolf2=NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 1/dim(otodat[otodat$source == 'ralston and miyamoto', ])[1], wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat)[1],  wt.lf = 1/length(pseudolf$curr_month_year), wt.lf2 = 0)
boot_time =  (proc.time() -  timer8)[3] / 60 / 60
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 8 complete!"))

print('Booting Model 9')
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 9")
timer9 = proc.time()
bootstrap_results$booted_param_ests_model9 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_9', boot_iterations = boot_iterations, tagdat=tagdat,  otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], pseudolf=pseudo_data, pseudolf2=NULL, wt.oto= 1, wt.oto2= 1, wt.oto3=1, wt.oto4=1, wt.tag = 1,  wt.lf = 1, wt.lf2 = 0)
boot_time =  (proc.time() -  timer9)[3] / 60 / 60
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 9 complete!"))

print('Booting Model 10')
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 10")
timer10 = proc.time()
bootstrap_results$booted_param_ests_model10 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_10', boot_iterations = boot_iterations, tagdat=tagdat, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], pseudolf=pseudo_data, pseudolf2 = NULL, wt.oto= 1/dim(otodat[otodat$source == 'demartini', ])[1], wt.oto2= 0, wt.oto3=1/dim(otodat[otodat$source == 'andrews bomb carbon', ])[1], wt.oto4=1/dim(otodat[otodat$source == 'andrews lead radium', ])[1], wt.tag = 1/dim(tagdat)[1], wt.lf = 1/length(pseudolf$curr_month_year), wt.lf2 = 0)
boot_time =  (proc.time() -  timer10)[3] / 60 / 60
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 10 complete!"))

print('Booting Model 11')
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "model 11")
timer11 = proc.time()
bootstrap_results$booted_param_ests_model12 = bootstrap_growth_params(filename = 'bootstrapped_parameter_estimates_model_11', boot_iterations = boot_iterations, tagdat=tagdat, tagdat2 = tagdat2, otodat=otodat[otodat$source == 'demartini', ], otodat2=otodat[otodat$source == 'ralston and miyamoto', ], otodat3=otodat[otodat$source == 'andrews bomb carbon', ], otodat4=otodat[otodat$source == 'andrews lead radium', ], pseudolf=pseudo_data, pseudolf2=NULL, wt.oto= 1, wt.oto2= 0, wt.oto3=1, wt.oto4=1, wt.tag = 1, wt.tag2 = 1, wt.lf = 1, wt.lf2 = 0)
boot_time =  (proc.time() -  timer11)[3] / 60 / 60
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste(round(boot_time, digits = 2), "Hours later, bootstrapping model 11 complete!"))
