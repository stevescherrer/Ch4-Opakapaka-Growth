#### Length Frequency Breakdown

## Two age classes
## Spawn in September t0 = september
## Start analysis in october
## Recruit at 2 months
## Followed for 17 Months
## SD up to 25% of individual's length

## Need to end up in the format:
# mode.age, mode.length, se.mode

install.packages('mixtools')
library('mixtools')
library('lubridate')

von_b_eq = function(t, t0 = -0.33, linf = 65.955, k = 0.2369){
  ## Get estimated length at time t using von bertalanffy function
  return(linf*(1-exp(-k*(t-t0))))
}

t1 = seq(0.083333, 1.416667, by = 1/12)
t2 = seq(1.08333, 2, by = 1/12)

## Dummy catch data - Made up of age and length data
dcd = data.frame()

## Fish lengths
for(i in 1:length(t1)){
  fish_lengths_1 = rnorm(100, mean(von_b_eq(t1)), sd = (abs(von_b_eq(t1)) * .25))
  fish_1 = sample(fish_lengths_1, round(rnorm(30, 20)[1]))
  fish_lengths_2 = rnorm(100, mean(von_b_eq(t2)), sd = (abs(von_b_eq(t1)) * .25))
  fish_2 = sample(fish_lengths_2, round(rnorm(30, 20)[1]))
  dcd = rbind(dcd, rbind(cbind(i, fish_1, 1), cbind(i, fish_2, 2)))
}
colnames(dcd) = c('month', 'length', 'year')
hist(dcd$length[dcd$month == 1], breaks = seq(0, 40, 2))

months = c('oct','nov','dec', "jan","feb","mar",'apr','may','jun','jul','aug','sep', 'oct','nov','dec')
dates = c(paste(months[1:3], '1989'), paste(months[4:length(months)], '1990'), paste(months[4:5], '1991'))
dcd$month = dates[dcd$month]
dcd = dcd[-which(dcd$month %in% dates[c(3, 8, 10, 15)]), ]


lfdat = data.frame()
for(i in 1:length(unique(dcd$month))){
  month_data = dcd[dcd$month == unique(dcd$month)[i], ]
  mode.age = which(dates == unique(dcd$month)[i]) / 12
  mode.age = c(mode.age, mode.age + 1)
  if(min(mode.age) > 1){
    mode.age = mode.age - 1
  }
  decomp = normalmixEM(month_data$length, k = 2)
  mode.len = decomp$mu
  mode.se = c(decomp$lambda * dim(month_data)[1])
  est.n = c(decomp$lambda * dim(month_data)[1])
  lfdat = rbind(lfdat, cbind(mode.age, mode.len, mode.se))
}

plot(x = lfdat$mode.age, y = lfdat$mode.len)


## Example run
p0 <- c(  70,   5.0, .10,   1.0,   .10,    5.0,      0,  0,   10.0,   10.0)
lb <- c(   50,   1.0, .05,    .1,   .05,    1.0,      0, -2,    0.0,    0.0)
ub <- c(  100,  15.0, .40,   1.5,   .50,   10.0,      0,  1,   15.0,   15.0)
fit.vb <- nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=tagdat,otodat=otodat,lfdat=lfdat)
fit.vb




### Now decomposing each month

