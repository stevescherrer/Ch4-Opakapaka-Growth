#### Making Plots for Dissertation Defense

### Keynote color palette
black = '#222222'
blue = '#39A7D8'
grey = '#838788'
white = '#A6AAA9'
green = '#8BBE62'
red = '#E74A43'

library(plotrix)

###### Generating Images For Growth Paper
proj_dir = '/Volumes/GoogleDrive/My Drive/Weng Lab/Personal_Folders/Steve/dissertation work/Ch 4. Opakapaka Growth/Analysis'
results_dir = '/Volumes/GoogleDrive/My Drive/Weng Lab/Personal_Folders/Steve/dissertation work/Ch 4. Opakapaka Growth/Presentations/Dissertation Defense'
fig_dir = '/Volumes/GoogleDrive/My Drive/Weng Lab/Personal_Folders/Steve/dissertation work/Dissertation/Defense Figures'


lit_vbgc_params = read.csv(file.path(proj_dir, 'data/lit_vbgf_params.csv'), stringsAsFactors = FALSE)
lit_vbgc_params = lit_vbgc_params[!is.na(lit_vbgc_params$Linf), ]
colnames(lit_vbgc_params) = c('author', 'n', 'linf', 'k', 't0', 'region', 'method')


#### Some general Functions
std_error = function(x){
  #### Calculates standard error of set (x)
  sqrt(var(x)/length(x))
}

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


##### Plotting Figures ####

#### Histogram of length at tagging and recapture 
pdf(file.path(fig_dir, 'Figure 1: Tag recapture length histogram.pdf'), height = 645, width = 500, unit = 'pix')
par(mfrow = c(3, 1), bg = black, col.axis = grey, col = grey, col.lab = grey, fg = grey)
hist(paka_growth$Lm, col = blue, border = black, xlim = c(20, 80), ylim = c(0, 200), xlab = 'Reported Fork Length (cm)', ylab = "n Individuals", main = NULL)
hist(paka_growth$Lr, col = blue, border = black, xlim = c(20, 80), ylim = c(0, 200), xlab = 'Reported Fork Length (cm)', ylab = "n Individuals", main = NULL)
# legend('topright', c("Tagged", "Recaptured"), fill = c(white, blue))
hist(paka_growth$dt, col = blue, border = black, xlim = c(0, 12), ylim = c(0, 100), xlab = 'Time at Liberty (years)', ylab = "n Individuals", main = NULL, breaks = seq(0, 11, .5))
dev.off()


#### Plotting predicted vs. observed length at recapture 
#### For our Modesl

### Constructing data frame of predicted values and residuals under each model 
predicted_recapture_lengths = data.frame(NULL)
predicted_recapture_residuals = data.frame(NULL)

models_to_predict = lit_vbgc_params[lit_vbgc_params$author %in%  "Maximum Likelihood - Integrative Model", ]
models_to_predict$linf = as.numeric(models_to_predict$linf)
models_to_predict$k = as.numeric(models_to_predict$k)
models_to_predict$t0 = as.numeric(models_to_predict$t0)
models_to_predict$t0[is.na(models_to_predict$t0)] = 0

for(i in 1:nrow(models_to_predict)){
  predicted_recapture_lengths = rbind(predicted_recapture_lengths, predict_recapture_length(Lm = tagdat[ ,1], dt = tagdat[ ,4], linf = as.numeric(models_to_predict$linf[i]), k = as.numeric(models_to_predict$k[i])))
  predicted_recapture_residuals = rbind(predicted_recapture_residuals, predict_recapture_length(Lm = tagdat[ ,1], dt = tagdat[ ,4], linf = as.numeric(models_to_predict$linf[i]), k = as.numeric(models_to_predict$k[i])) - tagdat[ ,2])
}

rownames(predicted_recapture_lengths)   = c('Model 11')
rownames(predicted_recapture_residuals) = c('Model 11')

pdf(file.path(fig_dir, 'Figure 3 - Predicted vs. Observed LR with validation data.pdf'), height = 8.5, width = 8.5)
par(mfrow = c(1, 1), bg = black, col.axis = grey, col.main = white, col = grey, col.lab = grey, fg = grey)
for(i in 1:nrow(predicted_recapture_lengths)){
  model_id = row.names(predicted_recapture_lengths)[i]
  plot(y = predicted_recapture_lengths[i, ], x = tagdat[ ,2],
       xlab = 'Observed Recapture FL (cm)', xlim = c(15, 80), 
       ylab = 'Predicted Recapture FL (cm)', ylim = c(15, 80),
       main = paste(model_id,'\n Linf = ', models_to_predict$linf[i],', K = ', models_to_predict$k[i], sep = ""),
       col = blue,
       pch = 19)
  abline(1, 1, lty = 2, col = white)
  arrows(x0 = 60, y0 = 70, x1 = 65, y1 = 70, length = 0.1, angle = 30, code = 2, col = white)
  text(x = 45, y = 70, labels = "Line of 1:1 agreement", cex = .75, col = white)
  # model_var = sum((predicted_recapture_lengths[i, ] - tagdat_validate[ ,2])^2) / length(tagdat_validate[ ,2])
  # text(x = 60, y = 30, labels = paste("Predictive Variance:", round(model_var, digits = 3)), cex = .75)
}
dev.off()


slices <- c(2, 1)
pdf(file.path(fig_dir, 'thirds pie chart.pdf'))
par(mfrow = c(1, 1), bg = black, col.axis = grey, col.main = white, col = grey, col.lab = grey, fg = grey)
pie3D(slices,explode=0.2,
      main="",
      col = c(blue, white),
      border = black)
### Note: Add and format labels for species in manually 
dev.off()





#### Barchart for Oct 1990 - Jan 1991 from Moffitt & Parish
### I ended up with one more record than they do but pretty close! (1048 vs. 1047).  number next to each month is the total number of fish for that month estimated from the histograms. * means this number was estimated a second time and matched
oct_1990 = data.frame('date' = as.POSIXct('1990-10-01'), 'val' = c(0, 0, 0, 0, 1, 0, 0, 2, 6, 17, 17, 15, 7, 5, 5, 0, 0, 0, 0), len = 6:24) # 75
nov_1990 = data.frame('date' = as.POSIXct('1990-11-01'), 'val' = c(0, 0, 0, 0, 0, 1, 1, 2, 0, 3, 3, 8, 5, 2, 1, 0, 0, 0, 0), len = 6:24) # 26 * 
jan_1991 = data.frame('date' = as.POSIXct('1991-01-01'), 'val' = c(0, 0, 0, 0, 0, 0, 0, 3, 2, 0, 12, 32, 30, 8, 3, 0, 0, 0, 0), len = 6:24) # 90 * 

oct = c()
for(i in 1:length(oct_1990$len)){
  oct = c(oct, rep(oct_1990$len[i], oct_1990$val[i]))
}

nov = c()
for(i in 1:length(nov_1990$len)){
  nov = c(nov, rep(nov_1990$len[i], nov_1990$val[i]))
}

jan = c()
for(i in 1:length(jan_1991$len)){
  jan = c(jan, rep(jan_1991$len[i], jan_1991$val[i]))
}

dev.off()


png(file.path(fig_dir, 'oct_1990.png'), height = 220, width = 384, units = "px")
par(mfrow = c(1, 1), bg = black, col.axis = grey, col.main = white, col = grey, col.lab = grey, fg = grey)
barplot(oct_1990$val, col = blue, names.arg = oct_1990$len, border = blue, ylim = c(0, 30))
dev.off()

png(file.path(fig_dir, 'nov_1990.png'), height = 220, width = 384, units = "px")
par(mfrow = c(1, 1), bg = black, col.axis = grey, col.main = white, col = grey, col.lab = grey, fg = grey)
barplot(nov_1990$val, col = blue, names.arg = nov_1990$len, border = blue, ylim = c(0, 30))
dev.off()


png(file.path(fig_dir, 'jan_1990.png'), height = 220, width = 384, units = "px")
par(mfrow = c(1, 1), bg = black, col.axis = grey, col.main = white, col = grey, col.lab = grey, fg = grey)
barplot(jan_1991$val, col = blue, names.arg = jan_1991$len, border = blue, ylim = c(0, 30))
dev.off()




#### Paka size von B
paka_size = von_b_eq(seq(0,15, .1), t0 = .219, linf = 67.6, k = .37)
  
pdf(file.path(fig_dir, 'paka_von_b.pdf'), width = 4.4, height = 3.5)
par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1, 1), bg = black, col.axis = grey, col.main = white, col = grey, col.lab = grey, fg = grey)
plot(FL..cm. ~ Age..Years., data = paka_size, col = blue, type = 'l', lwd = 4, ylab = 'Fork Length (cm)', xlab = 'Age (Years)')
dev.off()

von_b_eq = function(t, t0, linf, k){
  ## Get estimated length at time t using von bertalanffy function
  return(data.frame('Age (Years)' = t, 'FL (cm)' = linf*(1-exp(-k*(t-t0)))))
}






##### Random Figure VPN under different mortality rates

popi = 1000
b = .3
m = .10
f1 = .10
f2 = .40


estimate_population = function(n_init = 1000, b, m, f, iter = 10){
  n = c(n_init, rep(0, iter))
  for(t in 2:(iter+1)){
    n[t] = max(round(n[t-1]*(b - (m + f)) + n[t-1]), 0)
  }
  return(n)
}

age_struct = function(n_init, m, f, iter = 10){
  n = c(n_init, rep(0, iter))
  for(t in 2:(iter+1)){
    n[t] = max(round(n[t-1] - n[t-1]*(m + f)), 0)
  }
  return(n)
}


t = 1:15
nc = c(n_init, rep(0, 15))
for(t in 1:15){
  nc[t+1] = (f/(f+m))*n[t]*(1-exp(-(m + f)))
}

ypr = function(m, f = .1, n_init = 1000, iter = 15){
  nz = rep(0, iter)
  nc = rep(0, iter)
  n = c(n_init, rep(0, iter))
  for(t in 1:iter){
    nz[t] = (n[t]*(1-exp(-(m+f))))
    nc[t] = ((f/(f+m))*n[t]*(1-exp(-(m + f))))
    n[t+1] = n[t] - nz[t]
  }
  
  von_b_eq = function(t, t0, linf, k){
    ## Get estimated length at time t using von bertalanffy function
    return(linf*(1-exp(-k*(t-t0))))
  }
  paka_size = von_b_eq(1:iter, t0 = .219, linf = 67.6, k = .37)
  paka_weight = 0.0000381465*paka_size^2.79567
  
  return(sum(nc * paka_weight) / n_init)
}

yeilds = data.frame('f' = seq(0, .5, .01), 'ypr' = NA, stringsAsFactors = F)

for(i in 1:length(yeilds$f)){
  yeilds$ypr[i] = ypr(m, yeilds$f[i])
}


pdf(file.path(fig_dir, 'YPR.pdf'), height = 3.5, width = 4.4)
par(mfrow = c(1, 1), bg = black, col.axis = grey, col.main = white, col = grey, col.lab = grey, fg = grey)
plot(ypr ~ f, data = yeilds, type = 'l', ylab = 'Yeild Per Recruit (lbs)', xlab = 'Fishing Mortality (F)', col = blue, lwd = 4)
dev.off()




pop = data.frame('years' = 0:15)
pop$n_1 = estimate_population(n_init = 1000, b = b, m = m, f = f1, iter = 15)
pop$n_2 = estimate_population(n_init = 1000, b = b, m = m, f = f2, iter = 15)
pop$c_1 = age_struct(n_init = 1000, m = m, f = f1, iter = 15)
pop$c_2 = age_struct(n_init = 1000, m = m, f = f2, iter = 15)

pdf(file.path(fig_dir, 'Population Dynamics Under F.pdf'), height = 3.5, width = 4.4)
par(mfrow = c(1, 1), bg = black, col.axis = grey, col.main = white, col = grey, col.lab = grey, fg = grey)
plot(n_1 ~ years, data = pop, type = 'l', col = blue, lwd = 4, ylim = c(0, max(pop$n_1)), ylab = 'Population Size', xlab = 'Years')
lines(n_2 ~ years, data = pop, type = 'l', col = blue, lwd = 4)
dev.off()

pdf(file.path(fig_dir, 'Cohort Dynamics Under F.pdf'), height = 3.5, width = 4.4)
par(mfrow = c(1, 1), bg = black, col.axis = grey, col.main = white, col = grey, col.lab = grey, fg = grey)
plot(c_1 ~ years, data = pop, type = 'l', col = blue, lwd = 4, ylim = c(0, max(pop$c_1)), ylab = 'Cohort Size', xlab = 'Years Following Recruitment')
lines(c_2 ~ years, data = pop, type = 'l', col = blue, lwd = 4)
dev.off()

  t = 1:15
  nc = c(n_init, rep(0, 15))
  for(t in 0:15){
  nc[t+1] = (f/(f+m))*n[t+1](1-exp(-(m + f)))
  }
  