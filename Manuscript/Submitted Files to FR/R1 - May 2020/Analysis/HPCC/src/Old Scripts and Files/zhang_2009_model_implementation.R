#### Zhang et al 2009 implementation of Fabens Method
  ## Note: Similar to Francis but with additional A parameter.

install.packages('notifyR')
library('notifyR') # send_push()

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

laslett <- function(Linf, K, A){
  Lr_pred <- Lm + (Linf - Lm) * (1 - exp(-K*(A+dt)))
  dl_pred = Lr_pred - Lm
  dl_obs = Lr - Lm
  sigma = sqrt(sum(((dl_pred - dl_obs)^2)))
  NLL = -sum(log((1 / (sigma * sqrt(2*pi))) * (exp(-((dl_obs - dl_pred)^2 / (2*sigma^2))))))
  return(NLL)
}
laslett_mle <- mle2(laslett, start = list(Linf = l_inf_init, K = k_init, A = a_init), lower = list(Linf = 0.001, K = 0.001, A = 0.001), upper = list(Linf = 80, K = 1, A = 5), method = "L-BFGS-B")
summary(laslett_mle)

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
}

boot_param_ests = as.data.frame(boot_param_ests)
colnames(boot_param_ests) = c('Linf', 'K', 'A')
hist(boot_param_ests$A)


send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = "Boots Strapped!")





## How many fish fell within sigma of the data point
paka_growth$Lr_predicted = predict_recapture_length(Lm = paka_growth$Lm, dt = paka_growth$dt, linf = constant_mle@coef[1], k = constant_mle@coef[2])
paka_growth$dl_pred = paka_growth$Lr_predicted - paka_growth$Lm
paka_growth$sigma = sum(sqrt(((paka_growth$dl_pred - paka_growth$dL)^2)))/n

length(which(paka_growth$Lr >= paka_growth$Lr_predicted - paka_growth$sigma & paka_growth$Lr <= paka_growth$Lr_predicted + paka_growth$sigma)) / length(paka_growth$tag_id)
# 249 - 0.6484375

## Model AIC
AIC(constant_mle) # 3995.486


francis_constant_multiplier = function(Linf, K, nu, A){
  #print(c(Linf, K, nu))
  Lr_pred =  Lm + (Linf - Lm) * (1 - exp(-K*(A + dt)))
  dl_pred = Lr_pred - Lm
  dl_obs = Lr - Lm
  sigma = nu * dl_pred
  thing_to_exp = mpfr(-((dl_obs - dl_pred)^2 / (2*sigma^2)), precBits = 206)
  sig_rad = mpfr(sigma * sqrt(2*pi), precBits = 206)
  NLL = as.numeric(-sum(log((1 / sig_rad) * (exp(thing_to_exp)))))
  return(NLL)
}

## MLE2 WORKS HERE!
inverse_mle_2 <- mle2(francis_constant_multiplier,  start = list(Linf = l_inf_init, K = k_init, nu = nu_init,  A = a_init), lower = list(Linf = 0.001, K = 0.001, nu = 0.001, A = 0.001), method = "L-BFGS-B")
summary(inverse_mle_2)

# Coefficients:
#   Estimate Std. Error z value     Pr(z)    
# Linf 64.609955   2.301163 28.0771 < 2.2e-16 ***
#   K     0.197857   0.022096  8.9543 < 2.2e-16 ***
#   nu    0.421312   0.017731 23.7610 < 2.2e-16 ***
#   A     0.242328   0.042866  5.6532 1.575e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: 2028.955 

## How many fish fell within the sigma of the data point
paka_growth$Lr_predicted = predict_recapture_length(Lm = paka_growth$Lm, dt = paka_growth$dt, linf = inverse_mle@coef[1], k = inverse_mle@coef[2])
paka_growth$dl_pred = paka_growth$Lr_predicted - paka_growth$Lm
paka_growth$sigma = paka_growth$dl_pred * inverse_mle@coef[3]

length(which(paka_growth$Lr >= paka_growth$Lr_predicted - paka_growth$sigma & paka_growth$Lr <= paka_growth$Lr_predicted + paka_growth$sigma)) / length(paka_growth$tag_id)
# 239 - 0.6223958

### AIC
AIC(inverse_mle) # 2036.955


francis_exp_decline = function(Linf, K, nu, tau, A){
  # print(c(Linf, K, nu))
  Lr_pred =  Lm + (Linf - Lm) * (1 - exp(-K*(A +dt)))
  dl_pred = Lr_pred - Lm
  dl_obs = Lr - Lm
  sigma = tau * (1 - exp((-nu * dl_pred)))
  NLL = -sum(log((1 / (sigma * sqrt(2*pi))) * (exp(-((dl_obs - dl_pred)^2 / (2*sigma^2))))))
}

exp_decline_mle <- mle2(francis_exp_decline, start = list(Linf = 65.91, K = 0.23, nu = 0.09, tau = 4.3, A = a_init), lower = list(0.001, 0.001, 0.001, 0.001, 0.001), method = "L-BFGS-B")
summary(exp_decline_mle)

# Coefficients:
#   Estimate Std. Error z value     Pr(z)    
# Linf 65.960096   2.127064 31.0099 < 2.2e-16 ***
#   K     0.214789   0.021729  9.8850 < 2.2e-16 ***
#   nu    0.111559   0.021026  5.3056 1.123e-07 ***
#   tau   5.734874   0.646609  8.8692 < 2.2e-16 ***
#   A     0.124551   0.036333  3.4281 0.0006079 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: 1987.992 

paka_growth$Lr_predicted = predict_recapture_length(Lm = paka_growth$Lm, dt = paka_growth$dt, linf = exp_decline_mle@coef[1], k = exp_decline_mle@coef[2])
paka_growth$dl_pred = paka_growth$Lr_predicted - paka_growth$Lm
paka_growth$sigma = exp_decline_mle@coef[4] * (1 - exp(-exp_decline_mle@coef[3]  * paka_growth$dl_pred))

length(which(paka_growth$Lr >= paka_growth$Lr_predicted - paka_growth$sigma & paka_growth$Lr <= paka_growth$Lr_predicted + paka_growth$sigma)) / length(paka_growth$tag_id)
# 265 - 0.6901042

## AIC 
AIC(exp_decline_mle) 
# 1987.992 


francis_power = function(Linf, K, nu, tau, A){
  print(c(Linf, K, nu, A))
  Lr_pred =  Lm + (Linf - Lm) * (1 - exp(-K*(A + dt)))
  dl_pred = Lr_pred - Lm
  dl_obs = Lr - Lm
  sigma = nu*(dl_pred^tau)
  NLL = -sum(log(1 / (sqrt(2*pi) * sigma) * exp(-((dl_obs - dl_pred)^2 / (2*sigma^2)))))
  return(NLL)
}

power_decline_mle <- mle2(francis_power, start = list(Linf = 64.58, K = 0.252, nu = 1.113592, tau = 0.5060, A = a_init), lower = list(0.001, 0.001, 0.001, 0.001, 0.001), method = "L-BFGS-B")
summary(power_decline_mle)

# Coefficients:
#   Estimate Std. Error z value   Pr(z)    
# Linf 66.205675   2.097986 31.5568 < 2e-16 ***
#   K     0.220601   0.022107  9.9787 < 2e-16 ***
#   nu    0.970033   0.114992  8.4357 < 2e-16 ***
#   tau   0.562608   0.052766 10.6623 < 2e-16 ***
#   A     0.081289   0.037849  2.1477 0.03174 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: 1972.153 


paka_growth$Lr_predicted = predict_recapture_length(Lm = paka_growth$Lm, dt = paka_growth$dt, linf = power_decline_mle@coef[1], k = power_decline_mle@coef[2])
paka_growth$dl_pred = paka_growth$Lr_predicted - paka_growth$Lm
paka_growth$sigma = power_decline_mle@coef[3]*(paka_growth$dl_pred^power_decline_mle@coef[4])

length(which(paka_growth$Lr >= paka_growth$Lr_predicted - paka_growth$sigma & paka_growth$Lr <= paka_growth$Lr_predicted + paka_growth$sigma)) / length(paka_growth$tag_id)
# 266 - 0.6927083

## AIC
AIC(power_decline_mle) #  1972.153 

anova(constant_mle, inverse_mle_2, exp_decline_mle, power_decline_mle)













#### Laslett

## Notes: A = t1 - t0 = time at tagging
##        K = Growth coefficient


l1 = paka_growth$Lm
l2 = paka_growth$Lr
t1 = as.numeric(paka_growth$tm) / 60 / 60 / 24 / 365 # Convert from seconds to years
t2 = as.numeric(paka_growth$tr) / 60 / 60 / 24 / 365 # Convert from seconds to years

### Some starting values for debugging
Linf = 67.317556
K    = 0.213694
A    = 0.084914

laslett_likelihood = function(A, K, Linf){
  
  
  ### Length at capture
  ### Length at recapture
  ### t1
  ### t2
  
  
  ### Need to get out
    ## t
    ## A
    ## K
    ## Linf
  
  l1 = Linf * (1-exp(-K(A + t - t1)))
  l1
  
  
  
  ## Notes: Linf and A are assumed to be random variables. Linf assumed to have a normal distribution
  
  f = function(t1, t){
    dt = t - t1
    return(1 - exp(-K * (A + dt)))
  }
  
  ## Eq 1:
  l = function(t1, t){
    return(linf * f(t1, t))
  }

  t_from_l = function(l, t1){
    ## Function to return time t from length l
    t = (log((linf - l)/linf) / K) - A + t1
    return(t)
  }

  ### We also assume that measurements are taken at times t1 and t2 (with t1 < t2), so that
  ## Eq 2:  
    e1 = l1 - l(t1 = t1, t = t1)
    
  ## Eq 3:
    e2 = l2 - l(t1 = t2, t = t1)
    
    # where e1 and e2 represent measurement error. They are also normally distributed with mean 0 and variance σ2 and are in-
    # dependent from fish to fish. We assume also that e1 and e2 are independent of Linf and A.
    
    # We are now in a position to derive the joint distribution of
    # l1 and l2. For the moment, we argue conditional upon a known value of A, say A = a. Then l1 and l2 are both the sum of normal random variables and hence are themselves nor-
    # mal. Their first and second moments are µµ    
    mu1 = mean(l(t1, t) - l1)
    mu2 = mean(l(t2, t) - l2)
    
    sigma1 = var(l(t1, t) - l1)
    sigma2 = var(e2)
    
    cov_l1_l2 = cov(f)
  
  
  
  
  
  
  
  ## Linf and A are random and independent. Also that L inf is normally distributed
  
  dt = as.numeric(t2 - t1)
  
## Eq 1
## Equation 1b
f_t = function(t1){
  f = 1 - exp(-K * (A + t - t1))
  f[t1 - A > t] = 0
  return(f)
}


## Equation 1c
f = function(t1, t2){
  dt = t2 - t1
  return(1 - exp(-K * (A + dt)))
}

## Equation 1a
l = function(t){
  linf * f(t)
}


l1 = l(t1)
l2 = l(t2)

## Equations 2 and 3
e1 = l1 - l_t1
e2 = l2 - l_t2

mean1a = mean(f1(A, L1))
mean2a = mean(f2(A, L2))

sigma1a = var(f1(A, L1))
sigma2a = var(f2(A, L2))

covl1l2a = cov(f1(A, L1), f2(A, L2))
p_a = covl1l2a / (sigma1a * sigma2a)
q12a = ((l1 - mean1a)^2 / sigma1a^2) - 2*p_a * (((l1 - mean1a) * l2 - mean2a) / (sigma1a * sigma2a)) + ((l2-mean2a)^2 / sigma2a^2)

## Equation 4
hl1l2a = function(A){
  (1 / (2* pi * sigma1a * sigma2a * sqrt(1 - p_a^2))) * exp(-q12a / (2 * (1-p_a^2)))
}
## Equation 5
hl1l2 = integrateR(hl1l2a, lower = 0, upper = 1000000)
return(hl1l2)
  }
  




pg = paka_growth
paka_growth = paka_growth[1, ]

### Known values
## t1 = tagging time # in years
  t1 = as.numeric(paka_growth$tm) / 60 / 60 / 24 / 365 # in Years
## t2 = recapture time
  t2 = as.numeric(paka_growth$tr) / 60 / 60 / 24 / 365 # in Years
## dt = time increment between tagging and recapture
  dt = t2 - t1
## l1 = length at tagging
  l1 = paka_growth$Lm
## l2 = length at recapture
  l2 = paka_growth$Lr
## theta = a vector of fixed unknown parameters: = K
  theta = K

    ## t0 = theoretical age when fish is length 0
    ## A  = t1 - t0: varries from fish to fish. Age at tagging relative to t0 with density p(*)
    ## t  = real time
    ## linf = asymptotic length. Random from fish to fish with normal distribution
    ## mu_inf = mean of linf distribution
    ## sigma_inf = sigma of linf distribution

f = function(dt){
  1 - exp(-K * (A + dt))
}

l(t1)  
  
  
  
  
  ### Probability density function for A
  
  
  
  
    
    
    
t0 = 0
A = t1 - t0
# assumed growth is
l = function(t){
  linf * (1 - exp(-K * (A + t - t1)))
}

e1 = l(t1) - l1
e2 = l(t2) - l2

mu_1 = f(t1, t1)
mu_2 = f(t1, t2)
sigma1 = var(f(t1, t1))
sigma2 = var(ft(t1, t2))

covl112a = cov(f(t1, t1), f(t1, t2))
p_a = cor(l(t1), l(t2))

