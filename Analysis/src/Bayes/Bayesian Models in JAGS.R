### Running bayesian analysis with JAGS

proj_dir = getwd()
data_dir = file.path(proj_dir, "data")
src_dir = file.path(proj_dir, 'src')
results_dir = file.path(proj_dir, 'results')

# install.packages("R2jags")
library('R2jags')

### Loading in Datafile
mark_recapture_data = read.csv(file.path(data_dir, 'HO Mstr, temp (version 1).csv'), stringsAsFactors =  FALSE)
print(colnames(mark_recapture_data))
head(mark_recapture_data)

### Renaming data columns
colnames(mark_recapture_data) = c('tag_date', 'location', 'station', 'depth_f', 'species', 'previously_tagged', 'tag_id','fork_length_in', 'remarks', 'recapture_1_date', 'recapture_1_location', 'recapture_1_station', 'recapture_1_depth_f', 'recapture_1_fork_length_in', 'weight_1_lbs', 'days_1_free', 'growth_1_in', 'distance_1_miles','retagged_1',
                                  'recapture_2_date', 'recapture_2_location', 'recapture_2_station', 'recapture_2_depth_f', 'recapture_2_fork_length_in', 'weight_2_lbs', 'days_2_free', 'growth_2_in', 'distance_2_miles', 'retagged_2',
                                  'recapture_3_date', 'recapture_3_location', 'recapture_3_station', 'recapture_3_depth_f', 'recapture_3_fork_length_in', 'weight_3_lbs', 'days_3_free', 'growth_3_in', 'distance_3_miles', 'retagged_3',
                                  'recapture_4_date', 'recapture_4_location', 'recapture_4_station', 'recapture_4_depth_f', 'recapture_4_fork_length_in', 'weight_4_lbs', 'days_4_free', 'growth_4_in', 'distance_4_miles', 'x_retagged')


## How many total fish do we have in the data set?
n_tagged_fish = dim(mark_recapture_data)[1] # 4245!

### Subsetting out only Opakapaka with tag IDs - That is, fish that were marked
mark_recapture_data = mark_recapture_data[mark_recapture_data$species == '1' & mark_recapture_data$tag_id != '', ]
dim(mark_recapture_data)[1] # This gets you to the previously published 4179 tagged paka number from Kobayashi, Okamoto, & Oishi . for some reason doesn't exclude fish marked 'died'
n_tagged_paka = dim(mark_recapture_data)[1]

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


#### Restructuring data for JAGS 
## We know the maximum number of times a single fish was recaptured from the structure of our dataset
max_recaptures = 4
## Preallocating data objects in memory
L = matrix(nrow = length(mark_recapture_data$tag_id), ncol = (max_recaptures + 1), data = NA) # Note that matrix L gets an extra column corrosponding to initial fork length at tagging
dt = matrix(nrow = length(mark_recapture_data$tag_id), ncol = max_recaptures, data = NA) # Note that because dt is the difference in time, therefore it is one column less than L

## Defining a minimum time that must elapse before a recapture occurs
min_time_at_lib = 60 # in days

## Now lets reformat our data
for(i in 1:length(mark_recapture_data$tag_id)){
  ## For each individual, we start at the intial length at tagging and work our way up through recaptures
  L[i,1] = as.numeric(mark_recapture_data$fork_length_cm[i])
  ## Recapture 1
  if(!is.na(mark_recapture_data$recapture_1_fork_length_cm[i])){
    ## Calculate differece in time and if its greater than 60 days
    dt1 = difftime(mark_recapture_data$recapture_1_date[i], mark_recapture_data$tag_date[i], units = 'days') / 365
    if(all(!is.na(dt1) & dt1 > min_time_at_lib/365)){
      # Updating vector n
      # n[i] = n[i] + 1
      # Getting recorded fork length
      L[i, 2] = as.numeric(mark_recapture_data$recapture_1_fork_length_cm[i])
      # Calculating time at liberty
      dt[i, 1] = as.numeric(dt1)
    }
  }
  
  ### Repeat the above for subsequent recaptures
  ## Recapture 2
  if(!is.na(mark_recapture_data$recapture_2_fork_length_cm[i])){
    dt2 = difftime(mark_recapture_data$recapture_2_date[i], mark_recapture_data$tag_date[i], units = 'days') / 36
    if(all(!is.na(dt2) & dt2 > min_time_at_lib/365)){
     # n[i] = n[i] + 1
      L[i, 3] = as.numeric(mark_recapture_data$recapture_2_fork_length_cm[i])
      dt[i, 2] = as.numeric(dt2)
    }
  }
  
  ## Recapture 3
  if(!is.na(mark_recapture_data$recapture_3_fork_length_cm[i])){
    dt3 = difftime(mark_recapture_data$recapture_3_date[i], mark_recapture_data$tag_date[i], units = 'days') / 365
    if(all(!is.na(dt3) & dt3 > min_time_at_lib/365)){
      n[i] = n[i] + 1
      L[i, 4] = as.numeric(mark_recapture_data$recapture_3_fork_length_cm[i])
      dt[i, 3] = as.numeric(dt3)
    }
  }
  
  if(!is.na(mark_recapture_data$recapture_4_fork_length_cm[i])){
    dt4 = difftime(mark_recapture_data$recapture_4_date[i], mark_recapture_data$tag_date[i], units = 'days') / 365
    if(all(!is.na(dt4) & dt4 > min_time_at_lib/365)){
     # n[i] = n[i] + 1
      L[i, 5] = as.numeric(mark_recapture_data$recapture_4_fork_length_cm[i])
      dt[i,4] = as.numeric(dt4)
    }
  }
}

### Remove any fish missing an initial length
rm_ind = which(is.na(L[ ,1]))
L = L[-rm_ind, ]; dt = dt[-rm_ind, ] #, n = n[-rm_ind]; 

### L and dt matricies need to be adjusted for fish recaptured more than once but for whom one of their recaptures was excluded for being less than our 60 day threshold
## Set up vector indexing which fish to remove after our loop.We will remove fish that have no recaptures after 60 days a liberty.
valid_recaps = rep(TRUE, dim(dt)[1])

## Loop through each fish
for(i in 1:dim(dt)[1]){
  # If the first column of NAs does not come after the last column of data or if the whole row isn't NAs
  while(!all(is.na(dt[i, ])) && (!max(which(!is.na(dt[i, ]))) < min(which(is.na(dt[i, ]))))){
    # We shift all our data after the NA over by one and try again 
    L[i, (1 + min(which(is.na(dt[i, ])))):ncol(dt)] =  c(L[i, (1 + (min(which(is.na(dt[i, ]))))+ 1):ncol(dt)], NA)
    dt[i, min(which(is.na(dt[i, ]))):ncol(dt)] =  c(dt[i, (min(which(is.na(dt[i, ]))) + 1):ncol(dt)], NA)
    
  }
  if(all(is.na(dt[i, ]))){
    valid_recaps[i] = FALSE
  }
}

## Removing fish without valid recaptures
dt = dt[valid_recaps, ]
L = L[valid_recaps, ]

## Redetermining the number of unique capture/recapture events for each fish
n = rep(0, length(dt[,1]))
for(i in 1:length(n)){
  n[i] = length(L[i, ]) - length(which(is.na(L[i, ])))
}

## Check that the number of recaptures (n) for each fish is equal to the number of dts for each fish
for(i in 1:length(n)){
  if(any(is.na(L[i, 1:n[i]]) | is.na(dt[i, n[i]-1]))){
    print(i)
  }
}

## Wrap all of our data in a list for handing to JAGS
data = list(n = n, L = L, dt = dt, N = length(n))


#### Defining initial variables
## 3 lists corrosponding to each chain that will be initialized
inits1 = list('Linf_mu' = 50, 'Linf_tau' = 0.001, 'Shape' = 15, 'rate' = 9, 'k_mu' = 0.3, 'k_tau' = 2000, 'tau' = 0.4)
inits2 = list('Linf_mu' = 60*2, 'Linf_tau' = 0.001, 'Shape' = 1, 'rate' = 21, 'k_mu' = 0.3*2, 'k_tau' = 0.001, 'tau' = 0.001)
inits3 = list('Linf_mu' = 60/2, 'Linf_tau' = 0.001, 'Shape' = 1, 'rate' = 21, 'k_mu' = 0.3/2, 'k_tau' = 0.001, 'tau' = 0.001)
## Creating a single list that contains the lists corrosponding to each chain's initial values
inits = list(inits1, inits2, inits3)

#### Running our models
model_t = jags(data, inits, 
               model.file = file.path(src_dir,'Bayes/ModelTest_Jags.txt'),
               # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               parameters.to.save =  c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               DIC = T, 
               n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50)

model_1 = jags(data, inits, 
               model.file = file.path(src_dir,'Bayes/Model1_Jags.txt'),
               # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               parameters.to.save =  c('Linf_mu', 'Linf_std', 'k_mu', 'k_std'),
               DIC = T, 
               n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50)

model_2 = jags(data, inits, 
               model.file = file.path(src_dir,'Bayes/Model2_Jags.txt'),
               # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               parameters.to.save =  c('Linf_mu', 'Linf_std', 'k_mu', 'k_std'),
               DIC = T, 
               n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50)

model_3 = jags(data, inits, 
               model.file = file.path(src_dir,'Bayes/Model3_Jags.txt'),
               # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               parameters.to.save =  c('Linf_mu', 'Linf_std', 'k_mu', 'k_std'),
               DIC = T, 
               n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50)

model_4 = jags(data, inits, 
               model.file = file.path(src_dir,'Bayes/Model4_Jags.txt'),
               # parameters = c('Linf_std', 'k_std', 'var', 'Linf_mu', 'Linf_tau', 'shape', 'rate', 'k_mu', 'k_tau', 'tau'), 
               parameters.to.save =  c('Linf_mu', 'Linf_std', 'k_mu', 'k_std'),
               DIC = T, 
               n.chains = 3, n.iter = 500000, n.burnin = 10000, n.thin = 50)
