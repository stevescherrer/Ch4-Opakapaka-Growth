##### Fitting VBGF for Pristipomoides filamentosus with Mark Recapture Data
#### Written by: Stephen Scherrer with some code modified from Erik Franklin (2017)
#### Written on: 12 July 2017
#### Contact: scherrer@hawaii.edu
#### All Wrongs Preserved

##### Description ##### 
  #### A Script for fitting von Bertalanffy Growth Curve Parameters to mark - recapture 
  #### data collected by Henry Okimoto and DAR between 1989 and 1993 (see Okamoto 1993 for methods)
  #### using the following fitting methods:
  ####   1. Faben's Method
  ####   2. Haddon's Method for Faben
  ####   3. Wang's Method # 1
  ####   4. Wang's Method # 2
  ####   5. MLE - Constant Variance
  ####   6. MLE - Inverse Linear Between sd and E(deltaL)

  #### Script also plots VBGF derived from parameter estimates from these 6 methods and prior literature

  #### Note: Data Contains 487 recaptures were reported of 431 Individuals tagged between "1989-10-25 HST" and "1994-11-01 HST"
  ####  With Recaptures dating to "2003-10-12 HST". Removal of negative growth and capture history with dt < 60 days results in 384 data points

##### Script Set Up #####
#### Clearing Worksapce
rm(list = ls())

#### Setting a Script Timer
script_timer = proc.time()

#### Declaring Directory Path
src_dir = file.path(getwd(), 'src')
data_dir = file.path(getwd(), 'data')
results_dir = file.path(getwd(), 'results')
figure_dir = file.path(getwd(), 'figure')

#### Installing Principle Dependencies
# install.packages('gdata')
library('gdata') # read.xls()
# install.packages('FSA')
library('FSA') # vbFuns()
# install.packages('nlstools')
library('nlstools') # 
# install.packages('notifyR')
library('notifyR') # send_push()
# install.packages('minpack.lm')
library('minpack.lm') # nlsLM
# install.packages('fishmethods')
library('fishmethods')
# install.packages('bbmle')
library('bbmle') # MLE()


bonefish_growth = read.csv('/Volumes/GoogleDrive/My Drive/Weng Lab/Personal_Folders/Steve/dissertation work/Ch 4. Opakapaka Growth/Analysis/data/bonefish_growth.csv')

    ### Constructing derived variable dl (change in growth)
    bonefish_growth$dL = bonefish_growth$Lr - bonefish_growth$Lm
    ### Removing any fish that have a NA value for dL or dt (There is a single fish which had no tagging length and 7 fish with no recapture dates)
    bonefish_growth = bonefish_growth[!is.na(bonefish_growth$dL) & !is.na(bonefish_growth$dt), ]
    
  #### Creating a subset data frame that removes recording errors in length and time
    bonefish_growth = subset(bonefish_growth, dL >= 0)
    bonefish_growth = subset(bonefish_growth, dt >= 60/365)
    
    bonefish_growth$aproximate_rate = bonefish_growth$dL/bonefish_growth$dt
    
    
##### Functions #####
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
    
##### Model Fitting #####
    params = data.frame(stringsAsFactors = FALSE)
 #### 1. Fitting a growth curve to tagging data using the traditional Faben method
  ### Selecting intial starting parameters
    l_inf_init = max(bonefish_growth$Lr)
    k_init = mean(bonefish_growth$dL / bonefish_growth$dt)
    t0_init = 0

    ### Fitting Faben's growth function to data using non-linear least squares method
    Fabens.sv = list(Linf = l_inf_init, K = k_init)
    fvb = vbFuns("Fabens")
    FVB1 = nls(Lr ~ fvb(Lm, dt, Linf, K), start = Fabens.sv, data = bonefish_growth)
    summary(FVB1)
    # Parameters:
    #   Estimate Std. Error t value Pr(>|t|)    
    # Linf 65.95540    1.47979   44.57   <2e-16 ***
    #   K     0.23691    0.01632   14.52   <2e-16 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 3.505 on 430 degrees of freedom
    params = rbind(params, data.frame('method' = 'Faben', 'K' = coefficients(FVB1)['K'], 'Linf' = coefficients(FVB1)['Linf']))
    
    bonefish_growth$Lr_predicted = predict_recapture_length(Lm = bonefish_growth$Lm, dt = bonefish_growth$dt, linf = coef(FVB1)[1], k = coef(FVB1)[2])
    bonefish_growth$dl_pred = bonefish_growth$Lr_predicted - bonefish_growth$Lm
    bonefish_growth$sigma = sum(sqrt((bonefish_growth$dl_pred - bonefish_growth$dL)^2))/n
    
    length(which(bonefish_growth$Lr >= bonefish_growth$Lr_predicted - bonefish_growth$sigma & bonefish_growth$Lr <= bonefish_growth$Lr_predicted + bonefish_growth$sigma))
    # 249
    

#### Gulland and Holt (1959) Method
# Calculating mean size of individuals
# bonefish_growth = bonefish_growth[bonefish_growth$aproximate_rate > 6, ]
mean_size = (bonefish_growth$Lr + bonefish_growth$Lm) / 2
growth_rates = ((bonefish_growth$Lr - bonefish_growth$Lm) / (bonefish_growth$dt))

mean_size = mean_size[-277]
growth_rates = growth_rates[-277]
gh_mod = lm(mean_size ~ growth_rates)
summary(gh_mod)
plot(mean_size ~ growth_rates, pch = 19, cex = 0.5, xlab = "growth rate (cm/year)", ylab = "mean size (cm)", main = "Gulland & Holt Growth Parameter Estimations", xlim = c(0, 100), ylim = c(0, 70))
points(mean_size[which(bonefish_growth$aproximate_rate <= 6)] ~ growth_rates[aproximate_rate <= 6], col = 'blue')
abline(gh_mod)
xint = -1 * (coefficients(gh_mod)[1]/coefficients(gh_mod)[2])
text(x = 40, y = max(mean_size) - 5, labels = c(paste("L infinity =", round(xint, digits = 3), "\n", 'K = ', -1*round(coefficients(gh_mod)[2], digits = 3))))



### leverage plot?
plot(gh_mod)




##### HADDON's Implementation of Francis 1988

#### Francis 1988 implementation of Fabens Method
  ## Using NLL approach from Haddon
# bonefish_growth = bonefish_growth[bonefish_growth$dL != 0,]
n = dim(bonefish_growth)[1]
Lm = bonefish_growth$Lm
Lr = bonefish_growth$Lr
dt = bonefish_growth$dt
dl = bonefish_growth$dL

l_inf_init = max(bonefish_growth$Lr)
k_init = 0.3
nu_init = 0.1
tau_init = 4.1

francis_constant_var <- function(Linf, K){
  Lr_pred <- Lm + (Linf - Lm) * (1 - exp(-K*dt))
  dl_pred = Lr_pred - Lm
  dl_obs = Lr - Lm
  sigma = sqrt(sum(((dl_pred - dl_obs)^2)))
  NLL = -sum(log((1 / (sigma * sqrt(2*pi))) * (exp(-((dl_obs - dl_pred)^2 / (2*sigma^2))))))
  return(NLL)
}
constant_variance_mle <- mle2(francis_constant_var, start = list(Linf = l_inf_init, K = k_init), lower = list(0.001, 0.001), method = "L-BFGS-B")
summary(constant_variance_mle)


# Coefficients:
#   Estimate Std. Error z value     Pr(z)    
# Linf 62.246760   1.269855 49.0188 < 2.2e-16 ***
#   K     0.328333   0.033487  9.8048 < 2.2e-16 ***
#   ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: 1302.686 


## How many fish fell within sigma of the data point
bonefish_growth$Lr_predicted = predict_recapture_length(Lm = bonefish_growth$Lm, dt = bonefish_growth$dt, linf = constant_variance_mle@coef[1], k = constant_variance_mle@coef[2])
bonefish_growth$dl_pred = bonefish_growth$Lr_predicted - bonefish_growth$Lm
bonefish_growth$sigma = sum(sqrt(((bonefish_growth$dl_pred - bonefish_growth$dL)^2)))/n

length(which(bonefish_growth$Lr >= bonefish_growth$Lr_predicted - bonefish_growth$sigma & bonefish_growth$Lr <= bonefish_growth$Lr_predicted + bonefish_growth$sigma)) / length(bonefish_growth$tag_id)
# 0.6708861

## Model AIC
AIC(constant_variance_mle) # 1306.686


francis_inverse_linear = function(Linf, K, nu){
  #print(c(Linf, K, nu))
  Lr_pred =  Lm + (Linf - Lm) * (1 - exp(-K*dt))
  dl_pred = Lr_pred - Lm
  dl_obs = Lr - Lm
  sigma = nu * dl_pred
  thing_to_exp = mpfr(-((dl_obs - dl_pred)^2 / (2*sigma^2)), precBits = 206)
  sig_rad = mpfr(sigma * sqrt(2*pi), precBits = 206)
  NLL = as.numeric(-sum(log((1 / sig_rad) * (exp(thing_to_exp)))))
  return(NLL)
  }

## For some reason, this method provides estimates and std. error
inverse_linear_mle <- mle2(francis_inverse_linear,  start = list(Linf = l_inf_init, K = k_init, nu = nu_init), lower = list(0.001, 0.001, 0.001), method = "L-BFGS-B")
summary(inverse_linear_mle)

# Coefficients:
#   Estimate Std. Error
# Linf 75.1789712 2.16488719
# K     0.1375460 0.01583591
# nu    0.7337556 0.05954709
# 
# -2 log L: 655.0702 

## How many fish fell within the sigma of the data point
bonefish_growth$Lr_predicted = predict_recapture_length(Lm = bonefish_growth$Lm, dt = bonefish_growth$dt, linf = inverse_linear_mle@coef[1], k = inverse_linear_mle@coef[2])
bonefish_growth$dl_pred = bonefish_growth$Lr_predicted - bonefish_growth$Lm
bonefish_growth$sigma = bonefish_growth$dl_pred * inverse_mle@coef[3]

length(which(bonefish_growth$Lr >= bonefish_growth$Lr_predicted - bonefish_growth$sigma & bonefish_growth$Lr <= bonefish_growth$Lr_predicted + bonefish_growth$sigma)) / length(bonefish_growth$tag_id)
#  0.6708861

### AIC
AIC(inverse_linear) # 661.0702


francis_exp_decline = function(Linf, K, nu, tau){
  # print(c(Linf, K, nu))
  Lr_pred =  Lm + (Linf - Lm) * (1 - exp(-K*dt))
  dl_pred = Lr_pred - Lm
  dl_obs = Lr - Lm
 sigma = tau * (1 - exp((-nu * dl_pred)))
  thing_to_exp = mpfr(-((dl_obs - dl_pred)^2 / (2*sigma^2)), precBits = 206)
  sig_rad = mpfr(sigma * sqrt(2*pi), precBits = 206)
  NLL = as.numeric(-sum(log((1 / sig_rad) * (exp(thing_to_exp)))))
  }

#exp_decline_mle <- mle2(francis_exp_decline, start = list(Linf = 70, K = 0.23, nu = 0.09, tau = 4.3), lower = list(0.001, 0.001, 0.001, 0.001), method = "L-BFGS-B")
exp_decline_mle <- mle2(francis_exp_decline, start = list(Linf = 70, K = 0.23, nu = 0.09, tau = 4.3), lower = list(0.001, 0.001, 0.001, 0.001), method = "L-BFGS-B")

summary(exp_decline_mle)

# Coefficients:
#   Estimate Std. Error
# Linf 71.9241739 1.21712697
# K     0.1738155 0.01597755
# nu    0.3183782 0.10068976
# tau   3.3277715 0.66016535
# 
# -2 log L: 639.589 

bonefish_growth$Lr_predicted = predict_recapture_length(Lm = bonefish_growth$Lm, dt = bonefish_growth$dt, linf = exp_decline_mle@coef[1], k = exp_decline_mle@coef[2])
bonefish_growth$dl_pred = bonefish_growth$Lr_predicted - bonefish_growth$Lm
bonefish_growth$sigma = exp_decline_mle@coef[4] * (1 - exp(-exp_decline_mle@coef[3]  * bonefish_growth$dl_pred))

length(which(bonefish_growth$Lr >= bonefish_growth$Lr_predicted - bonefish_growth$sigma & bonefish_growth$Lr <= bonefish_growth$Lr_predicted + bonefish_growth$sigma)) / length(bonefish_growth$tag_id)
#  0.6265823

## AIC 
AIC(exp_decline_mle) 
# 647.589

francis_power = function(Linf, K, nu, tau){
  print(c(Linf, K, nu, tau))
  Lr_pred =  Lm + (Linf - Lm) * (1 - exp(-K*dt))
  dl_pred = Lr_pred - Lm
  dl_obs = Lr - Lm
  sigma = sign(dl_pred) * nu*abs(dl_pred)^tau
  thing_to_exp = mpfr(-((dl_obs - dl_pred)^2 / (2*sigma^2)), precBits = 206)
  sig_rad = mpfr(sigma * sqrt(2*pi), precBits = 206)
  NLL = as.numeric(-sum(log((1 / sig_rad) * (exp(thing_to_exp)))))
  return(NLL)
  }

power_decline_mle <- mle2(francis_power, start = list(Linf = 70, K = 0.23, nu = 0.09, tau = 4.3), lower = list(0.001, 0.1, 0.001, 0.1), method = "L-BFGS-B")
summary(power_decline_mle)

# Coefficients:
#   Estimate Std. Error
# Linf 70.9799264 1.40039051
# K     0.1777754 0.01735285
# nu    1.0144288 0.11731901
# tau   0.5798157 0.09946210
# 
# -2 log L: 635.2407 


bonefish_growth$Lr_predicted = predict_recapture_length(Lm = bonefish_growth$Lm, dt = bonefish_growth$dt, linf = power_decline_mle@coef[1], k = power_decline_mle@coef[2])
bonefish_growth$dl_pred = bonefish_growth$Lr_predicted - bonefish_growth$Lm
bonefish_growth$sigma = power_decline_mle@coef[3]*(bonefish_growth$dl_pred^power_decline_mle@coef[4])

length(which(bonefish_growth$Lr >= bonefish_growth$Lr_predicted - bonefish_growth$sigma & bonefish_growth$Lr <= bonefish_growth$Lr_predicted + bonefish_growth$sigma)) / length(bonefish_growth$tag_id)
# 0.6582278

## AIC
AIC(power_decline_mle) # 643.2407


anova(constant_variance_mle, inverse_linear_mle, exp_decline_mle , power_decline_mle)

## By AIC, 
AIC(constant_variance_mle) # 1306.686
AIC(inverse_linear_mle) # 646.3151
AIC(exp_decline_mle) # 647.589
AIC(power_decline_mle) # 643.2407

# Likelihood Ratio Tests
# Model 1: constant_variance_mle, [francis_constant_var]:
#   Linf+K
# Model 2: inverse_linear_mle, [francis_inverse_linear]:
#   Linf+K+nu
# Model 3: exp_decline_mle, [francis_exp_decline]: Linf+
#   K+nu+tau
# Model 4: power_decline_mle, [francis_power]: Linf+K+
#   nu+tau
# Tot Df Deviance    Chisq Df Pr(>Chisq)    
#   1      2  1302.69                           
#   2      3   655.07 647.6157  1  < 2.2e-16 ***
#   3      4   639.59  15.4812  1  8.333e-05 ***
#   4      4   635.24   4.3483  0  < 2.2e-16 ***
#   ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#### So the best model is the power decline model (model 4). Previously we fit it without respect to sex, now we'll break itout by sex and compare using Kaimura's likelihood ratio test
bonefish_growth_males = bonefish_growth[bonefish_growth$Sex == 'Male', ]
n = dim(bonefish_growth_males)[1]
Lm = bonefish_growth_males$Lm
Lr = bonefish_growth_males$Lr
dt = bonefish_growth_males$dt
dl = bonefish_growth_males$dL

power_decline_mle_males <- mle2(francis_power, start =  list(Linf = 70, K = 0.23, nu = 0.09, tau = 4.3), lower = list(30, 0.1, 0.001, .1), method = "L-BFGS-B")
summary(power_decline_mle_males)


bonefish_growth_males = bonefish_growth[bonefish_growth$Sex == 'Female', ]
n = dim(bonefish_growth_males)[1]
Lm = bonefish_growth_males$Lm
Lr = bonefish_growth_males$Lr
dt = bonefish_growth_males$dt
dl = bonefish_growth_males$dL

power_decline_mle_females <- mle2(francis_power, start = list(Linf = 70, K = 0.23, nu = 0.09, tau = 4.3), lower = list(0.001, 0.1, 0.001, 0.001), method = "L-BFGS-B")
summary(power_decline_mle_females)
