#### Laslett et al 2002 implementation of Fabens Method
## Note: Similar to Francis but with additional A parameter.

# install.packages('gdata')
library('gdata') # read.xls()
# install.packages('bbmle')
library('bbmle') # MLE()
# install.packages('notifyR')
library('notifyR') # send_push()
# install.packages('Deriv')
library('Deriv')
# install.packages('Rmpfr')
library('Rmpfr')
# install.packages('doParallel')
library('doParallel')
# install.packages('beepr')
library('beepr')
# install.packages('notifyR')
library('notifyR')
# install.packages('Rmpfr')
library('Rmpfr')

registerDoParallel(cores = detectCores()-1)

### Path to Data File
data_dir = 'data'

##### Loading and Cleaning Data Files #####
#### Mark Recapture Data
## Note: This data set is a fucking mess. Each line is a fish with location of capture, tag date, tagging depth, species, tag id, fork length, remarks, and then duplicate data for all of that each time it was recaptured up to 4 recaptures
mark_recapture_data = read.xls(file.path(data_dir, 'HO Mstr, temp (version 1).xlsx'), stringsAsFactors =  FALSE)

## How many total fish do we have in the data set?
dim(mark_recapture_data)[1] # 4245!

### Renaming data columns
colnames(mark_recapture_data) = c('tag_date', 'location', 'station', 'depth_f', 'species', 'previously_tagged', 'tag_id','fork_length_in', 'remarks', 'recapture_1_date', 'recapture_1_location', 'recapture_1_station', 'recapture_1_depth_f', 'recapture_1_fork_length_in', 'weight_1_lbs', 'days_1_free', 'growth_1_in', 'distance_1_miles','retagged_1',
                                  'recapture_2_date', 'recapture_2_location', 'recapture_2_station', 'recapture_2_depth_f', 'recapture_2_fork_length_in', 'weight_2_lbs', 'days_2_free', 'growth_2_in', 'distance_2_miles', 'retagged_2',
                                  'recapture_3_date', 'recapture_3_location', 'recapture_3_station', 'recapture_3_depth_f', 'recapture_3_fork_length_in', 'weight_3_lbs', 'days_3_free', 'growth_3_in', 'distance_3_miles', 'retagged_3',
                                  'recapture_4_date', 'recapture_4_location', 'recapture_4_station', 'recapture_4_depth_f', 'recapture_4_fork_length_in', 'weight_4_lbs', 'days_4_free', 'growth_4_in', 'distance_4_miles', 'x_retagged')

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

### Subsetting out only Opakapaka with tag IDs - That is, fish that were marked
mark_recapture_data = mark_recapture_data[mark_recapture_data$species == '1' & mark_recapture_data$tag_id != '', ]
dim(mark_recapture_data)[1] # This gets you to the previously published 4179 tagged paka number from Kobayashi, Okamoto, & Oishi . for some reason doesn't exclude fish marked 'died'

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
paka_growth = paka_growth[!is.na(paka_growth$dL) & !is.na(paka_growth$dt), ]

#### Creating a subset data frame that removes recording errors in length and time
paka_growth = subset(paka_growth, dL > 0)
paka_growth = subset(paka_growth, dt >= 60/365)


#### Modified Francis Approach ####
## Initializing data objects
# paka_growth = paka_growth[paka_growth$dL != 0,]
n = dim(paka_growth)[1]
Lm = paka_growth$Lm
Lr = paka_growth$Lr
dt = paka_growth$dt
dl = paka_growth$dL

l_inf_init = max(paka_growth$Lr)
k_init = mean(paka_growth$dL / paka_growth$dt)
nu_init = 0.1
tau_init = 4.1
a_init = 0

laslett_likelihood1 <- function(Linf, K, A){
  Lr_pred <- Lm + (Linf - Lm) * (1 - exp(-K*(A+dt)))
  dl_pred = Lr_pred - Lm
  dl_obs = Lr - Lm
  sigma = sqrt(sum(((dl_pred - dl_obs)^2)))
  NLL = -sum(log((1 / (sigma * sqrt(2*pi))) * (exp(-((dl_obs - dl_pred)^2 / (2*sigma^2))))))
  return(NLL)
}
laslett_mle1 <- mle2(laslett_likelihood1, start = list(Linf = l_inf_init, K = k_init, A = a_init), lower = list(Linf = 0.001, K = 0.001, A = 0.001), upper = list(Linf = 80, K = 1, A = 5), method = "L-BFGS-B")
summary(laslett_mle1)

# Coefficients:
#   Estimate Std. Error z value  Pr(z)    
# Linf 67.317556   2.127588 31.6403 <2e-16 ***
#   K     0.213694   0.024051  8.8852 <2e-16 ***
#   A     0.084914   0.065292  1.3005 0.1934    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: 3989.486 


#### Homemade bootstrapping since that isnt supported for MLE2 evidently
boot_intervals = 10000
boot_param_ests = as.data.frame(matrix(0, nrow = 0, ncol = 3))
# colnames(boot_param_ests) = c('Linf', 'K', 'A')

## Bootstrap sampling with replacement
# boot_param_ests = foreach(i = 1:boot_intervals, .combine = rbind) %dopar%{
for(i in 1:boot_intervals){
  ## Resampling data
  boot_sample = paka_growth[sample(nrow(paka_growth), size = nrow(paka_growth), replace = TRUE), ]
  
  ## Initializing data objects
  n = dim(boot_sample)[1]
  Lm = boot_sample$Lm
  Lr = boot_sample$Lr
  dt = boot_sample$dt
  dl = boot_sample$dL
  
  ## Refitting MLE Object
  boot_mle <- mle2(laslett, start = list(Linf = l_inf_init, K = k_init, A = a_init), lower = list(Linf = 0.001, K = 0.001, A = 0.001), upper = list(Linf = 80, K = 1, A = 10), method = "L-BFGS-B")
  
  ## Writing out parameters
  boot_param_ests = rbind(boot_param_ests, coef(boot_mle))
  # return(coef(boot_mle))
}

boot_param_ests = as.data.frame(boot_param_ests)
colnames(boot_param_ests) = c('Linf', 'K', 'A')
hist(boot_param_ests$Linf)
hist(boot_param_ests$K)
hist(boot_param_ests$A)

## Model AIC
AIC(constant_mle) # 3995.486



#### Laslett Likelihood Approach ####
Lr = paka_growth$Lr
Lm = paka_growth$Lm
t1 = paka_growth$tm
t2 = paka_growth$tr
dt = (t2 - t1) / 365

laslett_likelihood_5 = function(mu_inf, sigma_inf, K, mu_a){
  ## Estimating linf and a from length - Estimating the time tm - t0 of each individual based on its length
  pred_linf = function(){
    ## Note: Linfinity is normally distributed with a mean of mu_linf and variance sigma_sq_inf
    mpfr(Lr + (exp(-K * as.numeric(t2 - t1)/365) / (1-exp(-K * as.numeric(t2 - t1)/365))) * (Lr - Lm), precBits = 206)
  }
  
  pred_a = function(){
    ## Rearraning vonB to get time based on a length and a set of growth parameters
    mpfr((-1/K) * log(1-(Lm/pred_linf())), precBits = 206)
  }
  
  ### Let A = t1 – t0, then A varies from fish to fish
  # linf = pred_linf()
  # mu_inf = mpfr(mean(linf), precBits = 206)
  sigma_sq_inf = sigma_inf
  A = mu_a
  
  ## Some basic functions to describe the expected growth trajectory at times t1 and t2
  f1 = function(A){
    mpfr(1 - exp(-K * (A + as.numeric((t1 - t1))/365)), precBits = 206)
  }
  f2 = function(A){
    mpfr(1 - exp(-K * (A + as.numeric((t2 - t1))/365)), precBits = 206)
  }
  
  ## Getting the expected length of each individual at t11 and t2
  l1 = linf * f1(A)
  l2 = linf * f2(A)
  
  ## And the error associated
  e1 = l1 - Lm
  e2 = l2 - Lr
  sigma_sq = mpfr(sum(c(e1^2, e2^2))/length(c(e1, e2)), precBits = 206)
  
  a = A
  
  ### Calculating first and second moments of f1 and f2
  ## Mean around f1 and f2
  mu1 = function(a){
    mpfr(mu_inf * (f1(a)), precBits = 206)
  }
  mu2 = function(a){
    mpfr(mu_inf * (f2(a)), precBits = 206)
  }
  ## Variance around f1 and f2 
  sigma1_sq = function(a){
    mpfr((sigma_sq_inf * f1(a)^2 + sigma_sq), precBits = 206)
  }
  sigma2_sq = function(a){
    mpfr((sigma_sq_inf * f2(a)^2 + sigma_sq), precBits = 206)
  }
  
  ## Covariance of f1 f2
  cova = function(a){
    mpfr(sigma_sq_inf * f1(a) * f2(a), precBits = 206)
  }
  
  ## Some sort of correlation...
  rho = function(a){
   mpfr(cova(a) / (sqrt(sigma1_sq(a)) * sqrt(sigma2_sq(a))), precBits = 206)
  }
  
  ### Calculating the conditional density (in 2 parts)
  q12 = function(a){
    (((l1 - mu1(a))^2 / sigma1_sq(a)) - 2 * rho(a) * ((l1 - mu1(a)) * (l2 - mu2(a))) / (sqrt(sigma1_sq(a)) * sqrt(sigma2_sq(a))) + (l2 - mu2(a))^2 / sigma2_sq(a))
  }
  
  h = function(a){
    radical = mpfr(1 / (2 * pi * sqrt(sigma1_sq(a)) * sqrt(sigma2_sq(a)) * sqrt(1-rho(a)^2)), precBits = 206)
    thing_to_exp = mpfr(-q12(a) / (2*(1-rho(a)^2)), precBits = 206)
    return(as.numeric((radical) * exp(thing_to_exp)))
  }
  
  ## Calculating the unconditional joint density
  g = function(a){
    mpfr(h(a) * dlnorm(as.numeric(a), meanlog = as.numeric(mean(log(a))), sdlog = as.numeric(sd(log(a)))), precBits = 206)
  }
  
  
  g = function(a, i = NULL){
    if(is.null(i)){
      return(as.numeric(mpfr(h(a) * rho(a), precBits = 206)))
    }else{
      return(as.numeric(mpfr(h(a)[i] * rho(a)[i], precBits = 206)))
    }
  }
  
  
  # indv_joint_densities = c()
  # indv_joint_densities = foreach(i = 1:length(a), .combine = c) %dopar%{
  #   jd = integrate(f = Vectorize(g), lower = 0, upper = 1, i = 1)$value
  #   return(jd)
  # }
  
  NLL = -sum(log(indv_joint_densities))

   ## creating a function to optimize. The log of the g()
   log_g = function(a){
     mpfr(log(g(a)), precBits = 206)
   }
  
   ## Deriving likelihood
   ## Find max of log(g(a)) using brents method in optim (or whatever optimizer function)
   mu_g = function(a){
     ## Optimizing the function log(g(a))
     optimize(f = log_g, interval = c(as.numeric(mean(log(g(a))) - 1), as.numeric(mean(log(g(a))) + 1)), maximum = TRUE)$maximum
   }
  
   ## sigma_g: sesee formula in laslett. First fit a quadratic by choosing a value smaller than max, max, and larger than max, find funciton, double derive, then plug into formula
   sigma_g = function(a){
     x = c(mu_g(a) - 0.01, mu_g(a), mu_g(a) + 0.01)
     y = log(g(a))
     second_deriv_of_z = 2 * coef(lm(y ~ poly(x, 2, raw = TRUE)))[3]
     return(sqrt(-1 / second_deriv_of_z))
   }
  
     ## Then jump to equation 15. Note that last part of equation 15 is roughly sqrt(2*pi) - See second column.
   NLL = -sum(log(g(mu_g(a))) + log(sigma_g(a)) + log(sqrt(2*pi)))
  return(NLL)
}

### Some initial parameter values
mu_inf_init = 70
sigma_inf_init = 5
k_init = .1 # based on previous runs
mu_a_init = 1


## Estimating Model Parameters
laslett_mle5 <- mle2(laslett_likelihood_5, start = list(mu_inf = mu_inf_init, sigma_inf = sigma_inf_init, K = k_init, mu_a = mu_a_init), lower = list(mu_inf = 50, sigma_inf = 1, K = .5, mu_a = .1), upper = list(mu_inf = 110, sigma_inf = 15, K = .4, mu_a = 1.5), method = "L-BFGS-B")
summary(laslett_mle5)
beep(1)
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "Bootstrapping complete!")

##### Homemade bootstrapping since that isnt supported for MLE2 evidently.... #####
boot_mle = function(fun = laslett_likelihood4, data = paka_growth, start = list(mu_inf = l_inf_init, sigma_inf = .01, K = k_init, A = a_init), lower = list(mu_inf = 0.001, simga_inf = 0.001, K = 0.001, A = 0.001), upper = list(mu_inf = 100, sigma_inf = 0.001, K = 1, A = 10), method = "L-BFGS-B", boot_iterations = 10000){
  
  boot_param_ests = as.data.frame(matrix(0, nrow = 0, ncol = 4))
  # colnames(boot_param_ests) = c('Linf', 'K', 'A')
  
  ## Bootstrap sampling with replacement
  # boot_param_ests = foreach(i = 1:boot_intervals, .combine = rbind) %dopar%{
  for(i in 1:boot_intervals){
    # print(i)
    ## Resampling data
    boot_sample = data[sample(nrow(data), size = nrow(data), replace = TRUE), ]
    
    ## Initializing data objects
    
    ## Fish are tagged at time t1 with a length l1 and recaptured at time t2 with length l2
    t1 = boot_sample$tm[1]
    t2 = boot_sample$tr[1]
    Lm = boot_sample$Lm[1]
    Lr = boot_sample$Lr[1]
    
    ## Refitting MLE Object
    boot_mle <- mle2(fun, start = start, lower = lower, upper = upper, method = method)
    
    ## Writing out parameters
    boot_param_ests = rbind(boot_param_ests, coef(boot_mle))
    # return(coef(boot_mle))
  }
  boot_param_ests = as.data.frame(boot_param_ests)
  colnames(boot_param_ests) = c('Linf', 'K', 'A')
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "VBGF Script Complete! Run Time:seconds")
}

boot_param_ests = boot_mle(boot_iterations = 10)
hist(boot_param_ests$Linf)
hist(boot_param_ests$K)
hist(boot_param_ests$A)

#### Calling Home
beep(3)
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "Bootstrapping complete! (For realz this time)")

