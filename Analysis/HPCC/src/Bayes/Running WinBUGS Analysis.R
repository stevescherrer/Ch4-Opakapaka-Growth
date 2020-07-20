## Before running anything, need to first install WinBUGS on computer

## Workspace setup
setwd("G:/My Drive/Weng Lab/Data/Bottomfish/Okamoto_Mark_Recapture")
path2WinBugs = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/"
src_dir = file.path(getwd(), 'src')
results_dir = file.path(getwd(), 'results')
data_dir = file.path(getwd(), 'data')

## Install package dependencies 
# install.packages('R2WinBUGS')
library(R2WinBUGS)


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


## Read in BUGS Models

## Define data
L = cbind(tagdat[,1], tagdat[,2]) # L is a matrix with length = to the number of tagged fish where the first column is the marked length and the second is the length at recapture
dt = tagdat[,4] # t is the time at liberty in years
N = dim(L)[1]

mr_data = list('L', 'dt', 'N')

## Set initial parameters
inits1 = list(Linf_mu = 60, Linf_tau = 0.03, shape = 56, rate = 21, k_mu = 0.3, k_tau = 9121, tau = 0.12)#,
inits2 = list(Linf_mu = 60*.5, Linf_tau = 0.03, shape = 56, rate = 21, k_mu = 0.3 *0.5, k_tau = 9121, tau = 0.12)#,
inits3 = list(Linf_mu = 60*2, Linf_tau = 0.03, shape = 56, rate = 21, k_mu = 0.3 *2, k_tau = 9121, tau = 0.12)#,

inits = function(){
  inits1
}

## Run Models
model_1 = bugs(data = mr_data, inits = inits, 
               model.file = file.path(src_dir, 'Bayes/model1.bug.txt'),
              # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
              parameters.to.save =  c('Linf_mu', 'Linf_std', 'k_mu', 'k_std', 'Linf_tau', 'shape', 'rate', 'k_tau', 'tau'),
              DIC = T,
              n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))

model_2 = bugs(data, inits, 
               model.file = file.path(src_dir, 'model2.bug.txt'),
               parameters = c('Linf_std', 'k_std', 'var', 'k', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               n.chains = 1, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))

model_3 = bugs(data, inits, 
               model.file = file.path(src_dir, 'model3.bug.txt'),
               parameters = c('Linf_mu', 'Linf_std', 'Linf_tau', 'shape', 'k_mu', 'k_std', 'k_tau', 'rate', 'tau', 'var'), 
               n.chains = 1, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))

model_4 = bugs(data, inits, 
               model.file = file.path(src_dir, 'model4.bug.txt'),
               parameters = c('Linf_mu', 'Linf_std', 'Linf_tau', 'shape', 'k_mu', 'k_std', 'k_tau', 'rate', 'tau', 'var'), 
               n.chains = 1, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))


inits = function(){
  inits2
}

## Run Models
model_1.5 = bugs(data = mr_data, inits = inits, 
               model.file = file.path(src_dir, 'Bayes/model1.bug.txt'),
               # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               parameters.to.save =  c('Linf_mu', 'Linf_std', 'k_mu', 'k_std'),
               DIC = T,
               n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))

model_2.5 = bugs(data, inits, 
               model.file = file.path(src_dir, 'model2.bug.txt'),
               parameters = c('Linf_std', 'k_std', 'var', 'k', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               n.chains = 1, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))

model_3.5 = bugs(data, inits, 
               model.file = file.path(src_dir, 'model3.bug.txt'),
               parameters = c('Linf_mu', 'Linf_std', 'Linf_tau', 'shape', 'k_mu', 'k_std', 'k_tau', 'rate', 'tau', 'var'), 
               n.chains = 1, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))

model_4.5 = bugs(data, inits, 
               model.file = file.path(src_dir, 'model4.bug.txt'),
               parameters = c('Linf_mu', 'Linf_std', 'Linf_tau', 'shape', 'k_mu', 'k_std', 'k_tau', 'rate', 'tau', 'var'), 
               n.chains = 1, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))

inits = function(){
  inits3
}

## Run Models
model_1x2 = bugs(data = mr_data, inits = inits, 
               model.file = file.path(src_dir, 'Bayes/model1.bug.txt'),
               # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               parameters.to.save =  c('Linf_mu', 'Linf_std', 'k_mu', 'k_std'),
               DIC = T,
               n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))

model_2x2 = bugs(data, inits, 
               model.file = file.path(src_dir, 'model2.bug.txt'),
               parameters = c('Linf_std', 'k_std', 'var', 'k', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               n.chains = 1, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))

model_3x2 = bugs(data, inits, 
               model.file = file.path(src_dir, 'model3.bug.txt'),
               parameters = c('Linf_mu', 'Linf_std', 'Linf_tau', 'shape', 'k_mu', 'k_std', 'k_tau', 'rate', 'tau', 'var'), 
               n.chains = 1, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))

model_4x2 = bugs(data, inits, 
               model.file = file.path(src_dir, 'model4.bug.txt'),
               parameters = c('Linf_mu', 'Linf_std', 'Linf_tau', 'shape', 'k_mu', 'k_std', 'k_tau', 'rate', 'tau', 'var'), 
               n.chains = 1, n.iter = 500000, n.burnin = 10000, n.thin = 50,
               bugs.directory = "C:/Users/Weng Lab/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14",
               debug = TRUE, working.directory = file.path(results_dir, 'Bayes'))
