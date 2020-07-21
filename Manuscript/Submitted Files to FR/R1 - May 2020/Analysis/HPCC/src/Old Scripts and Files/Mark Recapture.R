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

  ### Subsetting out only Obonefishbonefish with tag IDs - That is, fish that were marked
  mark_recapture_data = mark_recapture_data[mark_recapture_data$species == '1' & mark_recapture_data$tag_id != '', ]
    dim(mark_recapture_data)[1] # This gets you to the previously published 4179 tagged bonefish number from Kobayashi, Okamoto, & Oishi . for some reason doesn't exclude fish marked 'died'

#### Now we want to format a table with lm (length at marking), lr (length at recapture), and dt (elapsed time)
 ### Note: If a fish was recaptured multiple times, there is a single entry for that individual corrosponding to the length at initial marking and the length at final recapture
    bonefish_growth = data.frame(stringsAsFactors = FALSE)
    for(i in 1:length(mark_recapture_data$tag_id)){
      if(!is.na(mark_recapture_data$recapture_4_fork_length_cm[i])){
        bonefish_growth = rbind(bonefish_growth, data.frame('tag_id' = mark_recapture_data$tag_id[i], 'Lm' = mark_recapture_data$fork_length_cm[i], 'Lr' = mark_recapture_data$recapture_4_fork_length_cm[i], 'tm' = mark_recapture_data$tag_date[i], 'tr' = mark_recapture_data$recapture_4_date[i], 'n_recaptures' = 4, stringsAsFactors = FALSE))
      }else if(!is.na(mark_recapture_data$recapture_3_fork_length_cm[i])){
        bonefish_growth = rbind(bonefish_growth, data.frame('tag_id' = mark_recapture_data$tag_id[i], 'Lm' = mark_recapture_data$fork_length_cm[i], 'Lr' = mark_recapture_data$recapture_3_fork_length_cm[i], 'tm' = mark_recapture_data$tag_date[i], 'tr' = mark_recapture_data$recapture_3_date[i], 'n_recaptures' = 3, stringsAsFactors = FALSE))
      }else if(!is.na(mark_recapture_data$recapture_2_fork_length_cm[i])){
        bonefish_growth = rbind(bonefish_growth, data.frame('tag_id' = mark_recapture_data$tag_id[i], 'Lm' = mark_recapture_data$fork_length_cm[i], 'Lr' = mark_recapture_data$recapture_2_fork_length_cm[i], 'tm' = mark_recapture_data$tag_date[i], 'tr' = mark_recapture_data$recapture_2_date[i], 'n_recaptures' = 2, stringsAsFactors = FALSE))
      }else if(!is.na(mark_recapture_data$recapture_1_fork_length_cm[i])){
        bonefish_growth = rbind(bonefish_growth, data.frame('tag_id' = mark_recapture_data$tag_id[i], 'Lm' = mark_recapture_data$fork_length_cm[i], 'Lr' = mark_recapture_data$recapture_1_fork_length_cm[i], 'tm' = mark_recapture_data$tag_date[i], 'tr' = mark_recapture_data$recapture_1_date[i], 'n_recaptures' = 1, stringsAsFactors = FALSE))
      }
    }
    bonefish_growth$dt = abs(difftime(bonefish_growth$tm, bonefish_growth$tr, units = "days"))    ## Converting dt from days to years
    bonefish_growth$dt = as.numeric(bonefish_growth$dt) / 365 # Converting to years
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
    
    
  #### 2. Model fitting by least squares - The way we did with Haddon back in MBIO 610
    ### This could totally be Francis... check it out..
    
    francis_mod = nls((dL ~ (l.inf - Lm) * (1-exp((-K * dt)))), data = bonefish_growth, 
                       start = list(K = k_init, l.inf = l_inf_init))
    summary(francis_mod)
    # Parameters:
    #   Estimate Std. Error t value Pr(>|t|)    
    #   K      0.23691    0.01632   14.52   <2e-16 ***
    #   l.inf 65.95546    1.47980   44.57   <2e-16 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 3.505 on 430 degrees of freedom
    params = rbind(params, data.frame('method' = 'Francis', 'K' = coefficients(francis_mod)['K'], 'Linf' = coefficients(francis_mod)['l.inf']))
    
    
  ### Model Diagnostics
    ssq = sum(resid(francis_mod)^2)
      # 5281.602
    plot(francis_mod)
    hist(resid(francis_mod))
    shapiro.test(resid(francis_mod))
      # Residuals appear evenly distributed around zero with increased variance for greater fitted values.  Shapiro test results beg to dissagree
    
    pdf(file.path(results_dir, 'francis_mod_residuals.pdf'))
    par(mfrow = c(2,1))
    ## Plot residuals vs length at release
    plot(nlsResiduals(francis_mod)$resi1[,2] ~ bonefish_growth$Lm, pch = 19, cex = .5, 
         main = "Residuals vs. FL at Release",
         xlab = 'Length at Release (cm)',
         ylab = 'Residuals')
    
    ## plot std residuals vs. length at release
    plot(nlsResiduals(francis_mod)$resi2[,2] ~ bonefish_growth$Lm, pch = 19, cex = .5, 
         main = "Standardized Residuals vs. FL at Release",
         xlab = 'Length at Release (cm)',
         ylab = 'Standardized Residuals')
    dev.off()
    
  #### 3. Fitting using Wang Method #1 (Wang 1998)
    wvb1 = vbFuns("Wang")
    Wang1.sv = list(Linf = l_inf_init, K = k_init, b = 0)
    WVB1 = nls(dL ~ wvb1(Lm, dt, Linf, K, b), start = Wang1.sv, data = bonefish_growth)
    summary(WVB1)
    
    # Parameters:
    #   Estimate Std. Error t value Pr(>|t|)    
    # Linf 66.83232    1.93130  34.605   <2e-16 ***
    #   K     0.22825    0.01935  11.797   <2e-16 ***
    #   b    -0.09405    0.12006  -0.783    0.434    
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 3.506 on 429 degrees of freedom
    # 
    params = rbind(params, data.frame('method' = 'Wang', 'K' = coefficients(WVB1)['K'], 'Linf' = coefficients(WVB1)['Linf']))

    
  #### 4. Fitting using Wang Method #2 (Wang 1998)
    wvb2 = vbFuns("Wang2")
    Wang2.sv = list(K = k_init, a = 200, d = 1)
    WVB2 = nls(dL ~ wvb2(Lm, dt, K, a, d), start = Wang2.sv, data = bonefish_growth)
    summary(WVB2)

    # Parameters:
    #   Estimate Std. Error t value Pr(>|t|)    
    #   K  0.22825    0.01935  11.797   <2e-16 ***
    #   a 69.91896    5.35203  13.064   <2e-16 ***
    #   d -1.09405    0.12006  -9.113   <2e-16 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 3.506 on 429 degrees of freedom
    params = rbind(params, data.frame('method' = 'Wang2', 'K' = coefficients(WVB2)['K'], 'Linf' = coefficients(WVB2)['a']))
    
    a = coefficients(WVB2)['a']
    d = coefficients(WVB2)['d']
    K = coefficients(WVB2)['K']
    X = bonefish_growth$Lm

    
###### The Rest of This Shit Was Shamelessly stolen from E. Franklin 2017 


  #### 5. EF: MLE constant variance
    ## Note: Constraining bonefish growth with aprox growth rates > 6 seems to fix fitting issues?
     # bonefish_growth = bonefish_growth[bonefish_growth$aproximate_rate > 6, ]

    n = dim(bonefish_growth)[1]
    Lm = bonefish_growth$Lm
    dt = bonefish_growth$dt
    dl = bonefish_growth$dL
    
    constant_var <- function(Linf, K){
      exp_delta_L <- (Linf - Lm) * (1 - exp(-K*dt))
      sigma <- sqrt(sum((exp_delta_L - dl)^2)/n)
      NLL <- -sum(dnorm(exp_delta_L, sigma, log = TRUE))
      return(NLL)
    }
    
    # francis_88 = function(Linf, K, s, m, p){
    # exp_delta_L <- (Linf - Lm) * (1 - exp(-K*dt))
    # sigma <- sqrt(sum((exp_delta_L - dl)^2)/n)
    # alpha_i = exp((-0.5 *(dl - exp_delta_L - m)^2 /  (sigma^2 + s^2) )/ (2*pi*(sigma^2 + s^2))^(1/2))
    # alpha = -sum(log((1-p)*alpha_i + p/(max(dl) - min(dl))))
    # return(alpha)
    # }
    # francis_mle = mle2(francis_88, start = list(Linf = 69, K = 0.22, s = .01, m = .01, p = .05))
    
    constant_mle <- mle2(constant_var, start = list(Linf = 69, K = 0.22))
    summary(constant_mle)
    # constant_mle@coef[1] -> Linf
    # constant_mle@coef[2] -> K
    # cons_fitted_values <- (Linf - Lm) * (1 - exp(-K*dt))

    
    nlm(f = constant_var, p = array(c(70,0.3), dim = c(2,1)))

    # Coefficients:
    #   Estimate Std. Error z value     Pr(z)    
    # Linf 46.943693   0.226379 207.367 < 2.2e-16 ***
    #   K     0.383439   0.012312  31.143 < 2.2e-16 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # -2 log L: 8032.019 

    

  #### 6. EF: MLE inverse linear between sd and E(deltaL)

    invline_var <- function(Linf, K, nu){
      exp_delta_L <- Lm + (Linf - Lm) * (1 - exp(-K*dt))
      sigma <- nu*exp_delta_L
      NLL <- -sum(dnorm(dl, exp_delta_L, sigma, log = TRUE))
      return(NLL)
    }
    
    invline_mle <- mle2(invline_var, start = list(Linf = 69, K = 0.22, nu = .3))
    
    
    summary(invline_mle)
    
    # Coefficients:
    #   Estimate Std. Error z value Pr(z)
    # Linf 53.94084         NA      NA    NA
    # K     0.86070         NA      NA    NA
    # nu    0.91898         NA      NA    NA
    # 
    # -2 log L: 3062.758 
    params = rbind(params, data.frame('method' = 'inverse_mle', 'K' = invline_mle@coef[2], 'Linf' = invline_mle@coef[1] ))
    
    
    invline_mle@coef[1] -> Linf
    invline_mle@coef[2] -> K
    invline_mle@coef[3] -> nu
    invl_fitted_values <- (Linf - Lm) * (1 - exp(-K*dt))
    
    bonefish_for_gtp = bonefish_growth
    colnames(bonefish_for_gtp) = c("tag_id", "L1", 'L2', 'T1', 'T2', 'n_recaptures', 'dt', 'dL')
    bonefish_for_gtp$T2 = as.numeric(abs(difftime(bonefish_for_gtp$T2, min(bonefish_for_gtp$T1), units = 'days')))
    bonefish_for_gtp$T1 = as.numeric(abs(difftime(bonefish_for_gtp$T1, min(bonefish_for_gtp$T1), units = 'days')))
    
    
    grotag(L1 = bonefish_for_gtp$L1, # Vector of length at release of tagged fish
           L2 = bonefish_for_gtp$L2, # Vector of length at recovery of tagged fish
           T1 = bonefish_for_gtp$T1, # Vector of julian time at release of tagged fish
           T2 = bonefish_for_gtp$T2, # Vector of julian time at recovery of tagged fish
           alpha = 35, # Numeric value giving an arbitrary length alpha
           beta = 50, # Numeric value giving an arbitrary length beta (beta > alpha)
           design = list(nu = 1, m = 0, p = 0, sea = 0), # List specifying the design of the model to estimate. Use 1 to designate whether a parameter(s) should be estimated. Type of parameters are: nu=growth vari- ability (1 parameter), m=bias parameter of measurement error (1 parameter), p=outlier probability (1 parameter), and sea=seasonal variation (2 parameters: u and w). Model 1 of Francis is the default settings of 0 for nu, m, p and sea.
           stvalue = list(sigma = 0.9, nu = 0.4, m = -1, p = 0.4, u = 0.4, w = 0.4), # Starting values of sigma (s) and depending on the design argument, nu, m, p, u, and w used as input in the nonlinear estimation (function optim) routine.
           upper = list(sigma = 5, nu = 1, m = 2, p = 0.5, u = 10, w = 10), # Upper limit of the model parameters’ (nu, m, p, u, and w) region to be investigated.
           lower = list(sigma = 0, nu = 0, m = -2, p = 0.1, u = 0, w = 0), # Lower limit of the model parameters’ (nu, m, p, u, and w) region to be investigated.
           gestimate = FALSE, # Logical specifying whether starting values of ga and gb (growth increments of alpha and beta) should be estimated automatically. Default = TRUE.
           st.ga = .6, # If gestimate=FALSE, user-specified starting value for ga.
           st.gb = .2, # If gestimate=FALSE, user-specified starting value for gb.
           st.galow = .1, # If gestimate=FALSE, user-specified lower limit for st.ga used in optimization.
           st.gaup = 1, # If gestimate=FALSE, user-specified upper limit for st.ga used in optimization.
           st.gblow = .1, # If gestimate=FALSE, user-specified lower limit for st.gb used in optimization.
           st.gbup = 1, # If gestimate=FALSE, user-specified upper limit for st.gb used in optimization.
           control = list(maxit = 100000000))
    
     grotagplus(tagdata = bonefish_for_gtp, 
               alpha = 0.1,
               beta = 0.2, 
               design = list(galpha=1,gbeta=1,s=1,nu=1,m=1,p=1,u=1,w=1),
               stvalue=list(s=0.81,nu=0.3,m=0,p=0.01,u=0.5,w=0.5),
               upper=list(s=3,nu=1,m=2,p=0.1,u=1,w=1),
               lower=list(s=0.1,nu=0.1,m=-2,p=0,u=0,w=0), debug = TRUE)

#### Model Selection ####
  alpha_threshold = 0.05
  if(summary(WVB1)$coefficients["b", 4] <= alpha_threshold){
    print("Wang's b is significant. Go with Wang!")
    vb_model = WVB1
  }else{
    print("Wang's b is insignficant. Go with Francis!")  
    vb_model = francis_mod
  }
    
    ## Comparing Wang's model to Francis's using ANOVA
    summary(WVB1)
    summary(francis_mod)
    summary(anova(WVB1, francis_mod))
    # Res.Df Res.Sum Sq Df  Sum Sq F value Pr(>F)
    # 1    429     5273.1                          
    # 2    430     5281.6 -1 -8.5015  0.6917 0.4061
    
    


#### Constructing Parameter Confidence Intervals through bootstrapping
    bootstrap_CI = nlsBoot(francis_mod, niter = 10000)
    
   
    


#### bonefish Summary Stats
## Number of fish tagged
length(unique(mark_recapture_data$tag_id))
  # 4172

## Number of fish recaptured
length(unique(bonefish_growth$tag_id))
  # 431

## percentage of recaptures
length(unique(bonefish_growth$tag_id)) / length(unique(mark_recapture_data$tag_id))
  # 0.1033078

## Number of each recaptures
  aggregate(bonefish_growth$tag_id, by = list(bonefish_growth$n_recaptures), FUN = length)
 # # of Recaps  # fish
 # 1 recap        394
 # 2 recap         35
 # 3 recap          1
 # 4 recap          2
  
## Range of sizes marked
fivenum(bonefish_growth$Lm)
  # 19.050 29.972 32.258 35.560 52.832

## Average length at marking
mean(bonefish_growth$Lm)
  # 32.8189

## SE of length at marking
se(bonefish_growth$Lm)
  # 0.2436386

hist(bonefish_growth$Lm, xlab = "Fork Length (cm)", main = "Length at Marking")


## Range of sizes recapture
fivenum(bonefish_growth$Lr)
  # 22.8600 36.0680 40.6400 46.8376 76.2000
mean(bonefish_growth$Lr)
  # 41.86461
se(bonefish_growth$Lr)
  # 0.4291505
hist(bonefish_growth$Lr, xlab = "Fork Length (cm)", main = "Length at Recapture")

## Time at liberty
## Days at liberty
fivenum(bonefish_growth$dt * 365)
  # 1.0  168.5  394.0  812.5 3748.0 ## Note: In Days
## Years at liberty
fivenum(bonefish_growth$dt)
  # 0.002739726  0.461643836  1.079452055  2.226027397 10.268493151
mean(bonefish_growth$dt * 365)
  # 600.125
se(bonefish_growth$dt * 365)
  # 29.96815

## Daily growth rate
fivenum(bonefish_growth$dL / (bonefish_growth$dt * 365))
  # -0.38100000  0.01080014  0.01657395  0.02327096  0.50800000
mean(bonefish_growth$dL / (bonefish_growth$dt * 365))
  # 0.01719882
se(bonefish_growth$dL / (bonefish_growth$dt * 365))
  # 0.001853376

## Time at Liberty vs. Growth
plot(bonefish_growth$dL ~ bonefish_growth$dt,
     main = "Growth vs. Time at Liberty",
     xlab = "Years",
     ylab = "dL (cm)")




## Negative Growth
mean(bonefish_growth$dL[bonefish_growth$dL < 0]) 
  # -0.8212667
se(c(bonefish_growth$dL[bonefish_growth$dL < 0], abs(bonefish_growth$dL[bonefish_growth$dL < 0])))^2
  # 0.03622536
mean(bonefish_growth$dt[bonefish_growth$dL < 0] * 365) 
  # 32.33333 days
se(bonefish_growth$dt[bonefish_growth$dL < 0] * 365) 
  # 5.604137 days

##### Diagnosing Measurement Error #####

### Method 1: Irrespective of dt
length(which(bonefish_growth$dL < 0))
  # 15
hist(bonefish_growth$dL[bonefish_growth$dL < 0], main = "Distribution of Negative Growth Measurements", xlab = 'cm')

fivenum(bonefish_growth$dL[bonefish_growth$dL < 0])
  # -2.540 -1.143 -0.762 -0.254 -0.254

measurement_error_variance_1 = sd(c(bonefish_growth$dL[bonefish_growth$dL < 0], abs(bonefish_growth$dL[bonefish_growth$dL < 0]))) ^2
  # 1.086761

### Plotting measurement error with time (dL ~ dT)
plot(dL ~ dt, data = bonefish_growth[bonefish_growth$dL < 0, ], pch = 19, main = "Measurement Error With Time")
lm.mod = lm(dL ~ dt, data = bonefish_growth[bonefish_growth$dL < 0, ])
abline(lm.mod)
dev.off()

### Method 2: With dt
neg_bonefish = bonefish_growth[bonefish_growth$dL < 0, ]
neg_bonefish$dt = neg_bonefish$dt / 365 # converting back to days from years
## Create daily estimates of bonefish growth size
growth_est = von_b_eq(t = 1:10000, t0 = 0, linf = params$Linf[params$method == 'Francis'], k = params$K[params$method == 'Francis']/365)
neg_growth = data.frame(stringsAsFactors = FALSE)
for(i in 1:length(neg_bonefish$dL)){
  tag_date_index = which.min(abs(growth_est - neg_bonefish$Lm[i]))
  tagging_length = growth_est[tag_date_index]
  est_recap_length = growth_est[tag_date_index + neg_bonefish$dt[i]]
  recap_length_index = which.min(abs(growth_est - neg_bonefish$Lr[i]))
  recap_length = growth_est[recap_length_index]
  est_tag_length = growth_est[recap_length_index - neg_bonefish$dt[i]]
  recap_error = neg_bonefish$dL[i] - (tagging_length - est_recap_length) 
  tagging_error =  neg_bonefish$dL[i] - (recap_length - est_tag_length)
  neg_growth = rbind(neg_growth, data.frame('tag_id' = neg_bonefish$tag_id[i], 'tagging_error' = tagging_error, 'recapture_error' = recap_error, stringsAsFactors = FALSE))
}

mean(c(neg_growth$tagging_error, neg_growth$recapture_error))
  # -0.8325608
measurement_error_variance_2 = sd(c(c(neg_growth$tagging_error, neg_growth$recapture_error), abs(c(neg_growth$tagging_error, neg_growth$recapture_error))))^2
  # 1.08804

##### End Of Script Clean Up #####
run_time = round((proc.time() - script_timer)[3])
send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste("VBGF Script Complete! Run Time:", run_time, "seconds"))





bonefish = bonefish_for_gtp
data(bonito)

#Model 4 of Francis (1988)
with(bonefish,
     grotag(L1=L1, L2=L2, T1=T1, T2=T2,alpha=35,beta=55,
            design=list(nu=1,m=1,p=1,sea=1),
            stvalue=list(sigma=0.9,nu=0.4,m=-1,p=0,u=0.4,w=0.4),
            upper=list(sigma=5,nu=1,m=2,p=0.5,u=1,w=1),
            lower=list(sigma=0,nu=0,m=-2,p=0,u=0,w=0),control=list(maxit=1e4)))




fun = vbFuns('Laslett')
fun_starts = 
WVB1 = nls(dL ~ wvb1(Lm, dt, Linf, K, b), start = Wang1.sv, data = bonefish_growth)
Wang1.sv = list(Linf = l_inf_init, K = k_init, b = 0)
summary(WVB1)


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


### Which are my outliers?
outliers = c(276, 36, 313, 275)
my_bad_fish = bonefish_growth[which(bonefish_growth$aproximate_rate <= 0.6), ]

c(which(bonefish_growth$dL < 1.086761), c(276, 36, 313, 275))








##### HADDON's Implementation

#### Francis 1988 implementation of Fabens Method
  ## Using NLL approach from Haddon
# bonefish_growth = bonefish_growth[bonefish_growth$dL != 0,]
n = dim(bonefish_growth)[1]
Lm = bonefish_growth$Lm
Lr = bonefish_growth$Lr
dt = bonefish_growth$dt
dl = bonefish_growth$dL

l_inf_init = max(bonefish_growth$Lr)
k_init = mean(bonefish_growth$dL / bonefish_growth$dt)
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
constant_mle <- mle2(francis_constant_var, start = list(Linf = l_inf_init, K = k_init), lower = list(0.001, 0.001), method = "L-BFGS-B")
summary(constant_mle)

# Coefficients:
#   Estimate Std. Error z value     Pr(z)    
# Linf 65.911119   1.639920  40.192 < 2.2e-16 ***
#   K     0.237526   0.018186  13.061 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: 3991.384 


## How many fish fell within sigma of the data point
bonefish_growth$Lr_predicted = predict_recapture_length(Lm = bonefish_growth$Lm, dt = bonefish_growth$dt, linf = constant_mle@coef[1], k = constant_mle@coef[2])
bonefish_growth$dl_pred = bonefish_growth$Lr_predicted - bonefish_growth$Lm
bonefish_growth$sigma = sum(sqrt(((bonefish_growth$dl_pred - bonefish_growth$dL)^2)))/n

length(which(bonefish_growth$Lr >= bonefish_growth$Lr_predicted - bonefish_growth$sigma & bonefish_growth$Lr <= bonefish_growth$Lr_predicted + bonefish_growth$sigma)) / length(bonefish_growth$tag_id)
# 249 - 0.6484375

## Model AIC
AIC(constant_mle) # 2069.652


francis_constant_multiplier = function(Linf, K, nu){
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
inverse_mle <- mle(francis_constant_multiplier,  start = list(Linf = l_inf_init, K = k_init, nu = nu_init), lower = list(0.001, 0.001, 0.001), method = "L-BFGS-B")
## This method uses MLE2 but can't provide standard error. Same estimates though.
inverse_mle_2 <- mle2(francis_constant_multiplier,  start = list(Linf = l_inf_init, K = k_init, nu = nu_init), lower = list(Linf = 0.001, K = 0.001, nu = 0.001), method = "L-BFGS-B")
summary(inverse_mle)

# Coefficients:
#   Estimate Std. Error
# Linf 58.5347590 1.01355047
# K     0.3461791 0.02024909
# nu    0.4596702 0.01966135
# 
# -2 log L: 2102.379 

## How many fish fell within the sigma of the data point
bonefish_growth$Lr_predicted = predict_recapture_length(Lm = bonefish_growth$Lm, dt = bonefish_growth$dt, linf = inverse_mle@coef[1], k = inverse_mle@coef[2])
bonefish_growth$dl_pred = bonefish_growth$Lr_predicted - bonefish_growth$Lm
bonefish_growth$sigma = bonefish_growth$dl_pred * inverse_mle@coef[3]

length(which(bonefish_growth$Lr >= bonefish_growth$Lr_predicted - bonefish_growth$sigma & bonefish_growth$Lr <= bonefish_growth$Lr_predicted + bonefish_growth$sigma)) / length(bonefish_growth$tag_id)
# 302 - 0.7864583

### AIC
AIC(inverse_mle) # 2076.66


francis_exp_decline = function(Linf, K, nu, tau){
  # print(c(Linf, K, nu))
  Lr_pred =  Lm + (Linf - Lm) * (1 - exp(-K*dt))
  dl_pred = Lr_pred - Lm
  dl_obs = Lr - Lm
  sigma = tau * (1 - exp((-nu * dl_pred)))
  NLL = -sum(log((1 / (sigma * sqrt(2*pi))) * (exp(-((dl_obs - dl_pred)^2 / (2*sigma^2))))))
}

exp_decline_mle <- mle2(francis_exp_decline, start = list(Linf = 65.91, K = 0.23, nu = 0.09, tau = 4.3), lower = list(0.001, 0.001, 0.001, 0.001), method = "L-BFGS-B")
summary(exp_decline_mle)

# Coefficients:
#   Estimate Std. Error z value     Pr(z)    
# Linf 64.010680   1.689689 37.8831 < 2.2e-16 ***
#   K     0.260592   0.020690 12.5954 < 2.2e-16 ***
#   nu    0.165102   0.024491  6.7414 1.569e-11 ***
#   tau   4.943619   0.366822 13.4769 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: 2006.085 

bonefish_growth$Lr_predicted = predict_recapture_length(Lm = bonefish_growth$Lm, dt = bonefish_growth$dt, linf = exp_decline_mle@coef[1], k = exp_decline_mle@coef[2])
bonefish_growth$dl_pred = bonefish_growth$Lr_predicted - bonefish_growth$Lm
bonefish_growth$sigma = exp_decline_mle@coef[4] * (1 - exp(-exp_decline_mle@coef[3]  * bonefish_growth$dl_pred))

length(which(bonefish_growth$Lr >= bonefish_growth$Lr_predicted - bonefish_growth$sigma & bonefish_growth$Lr <= bonefish_growth$Lr_predicted + bonefish_growth$sigma)) / length(bonefish_growth$tag_id)
# 285 - 0.7421875

## AIC 
AIC(exp_decline_mle) 
# 2014.085


francis_power = function(Linf, K, nu, tau){
  #print(c(Linf, K, nu))
  Lr_pred =  Lm + (Linf - Lm) * (1 - exp(-K*dt))
  dl_pred = Lr_pred - Lm
  dl_obs = Lr - Lm
  sigma = nu*(dl_pred^tau)
  NLL = -sum(log(1 / (sqrt(2*pi) * sigma) * exp(-((dl_obs - dl_pred)^2 / (2*sigma^2)))))
  return(NLL)
  }

power_decline_mle <- mle(francis_power, start = list(Linf = 65.91, K = 0.23, nu = 0.09, tau = 4.3), lower = list(0.001, 0.001, 0.001, 0.001), method = "L-BFGS-B")
summary(power_decline_mle)

# Coefficients:
#   Estimate Std. Error z value     Pr(z)    
# Linf 64.583930   1.678009  38.488 < 2.2e-16 ***
#   K     0.252336   0.018880  13.365 < 2.2e-16 ***
#   nu    1.113592   0.111052  10.028 < 2.2e-16 ***
#   tau   0.506026   0.044155  11.460 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: 1977.733 


bonefish_growth$Lr_predicted = predict_recapture_length(Lm = bonefish_growth$Lm, dt = bonefish_growth$dt, linf = power_decline_mle@coef[1], k = power_decline_mle@coef[2])
bonefish_growth$dl_pred = bonefish_growth$Lr_predicted - bonefish_growth$Lm
bonefish_growth$sigma = power_decline_mle@coef[3]*(bonefish_growth$dl_pred^power_decline_mle@coef[4])

length(which(bonefish_growth$Lr >= bonefish_growth$Lr_predicted - bonefish_growth$sigma & bonefish_growth$Lr <= bonefish_growth$Lr_predicted + bonefish_growth$sigma)) / length(bonefish_growth$tag_id)
# 278 - 0.78239583

## AIC
AIC(power_decline_mle) # 1985.733


anova(constant_mle, inverse_mle_2, exp_decline_mle, power_decline_mle)

## By AIC, 
AIC(constant_mle) # 2069.652
AIC(inverse_mle) # 2076.66
AIC(exp_decline_mle) # 1992.457
AIC(power_decline_mle) # 1985.733

# Likelihood Ratio Tests
# Model 1: constant_mle, [francis_constant_var]: Linf+K
# Model 2: inverse_mle, [francis_constant_multiplier]: Linf+K+nu
# Model 3: exp_decline_mle, [francis_exp_decline]: Linf+K+nu+tau
# Model 4: power_decline_mle, [francis_power]: Linf+K+nu+tau

# Tot Df Deviance  Chisq Df Pr(>Chisq)    
# 1      2  -121.44                         
# 2      3  -131.35 9.9140  1    0.00164 ** 
# 3      4  -132.58 1.2291  1    0.26759    
# 4      4  -132.72 0.1401  0    < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1





##### Can we apply shnute's model?

# yr = length recap
# ym = length mark

yr = bonefish_growth$Lr
ym = bonefish_growth$Lm
dt = bonefish_growth$dt
y2 = 72.1
tau_2 = 45.9

schnute = function(a, b, y1, tau1){
  yr_pred = (ym^b * exp(-a*(dt)) + (y2^b - y1^b * exp(-a * (tau2 - tau1))) * (1 - exp(-a * dt))/ 1-exp(-a * (tau2 - tau1)))^(1/b)
  sigma = sqrt(sum(((yr - yr_pred)^2)))
  NLL = -sum(log((1 / (sigma * sqrt(2*pi))) * (exp(-((yr_pred - yr)^2 / (2*sigma^2))))))
  return(NLL)
}


shnute_mle <- mle2(schnute, start = list(a = 0, b = 0, y1 = 1, tau1 = 3), lower = list(a = 0.01, b = 0.1, y1 = .01,  tau1 = 1), upper = list(a = 5, b = 5, y1 = 50,  tau1 = 5),  method = "L-BFGS-B")

a = coef(shnute_mle)[1]
b = coef(shnute_mle)[2]
y1 = coef(shnute_mle)[3]
tau1 = coef(shnute_mle)[4]

Lr_pred = (ym^b * exp(-a*(dt)) + (y2^b - y1^b * exp(-a * (tau2 - tau1))) * (1 - exp(-a * dt))/ 1-exp(-a * (tau2 - tau1)))^(1/b)

par(mfrow = c(2, 1))
plot(Lr_pred ~ bonefish_growth$Lr, ylim = c(0, 80), xlim = c(0,80))
abline(1,1)


plot((Lr_pred - bonefish_growth$Lr) ~ bonefish_growth$Lr)
abline(1,0)
