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
