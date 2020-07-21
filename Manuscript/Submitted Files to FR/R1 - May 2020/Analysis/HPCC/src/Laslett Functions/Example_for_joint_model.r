
# Simulate VB tag-recapture, otolith and length-frequency growth data
# =====================================================================

# Specify VB parameters common to all datasets(noting that Linf ~ Norm(mu.Linf, sigma.Linf) )
k <- 0.2
a0 <- 0
mu.Linf <- 100
sigma.Linf <- 5

# a)tag-recapture data:
# ----------------------
set.seed(1)
n <- 500

# generate log-normal release ages (relative to a0)
par1.A <- 0.5
par2.A <- 0.3
A <- rlnorm(n,par1.A,par2.A)

# generate times at liberty assuming gamma distribution
dt <- rgamma(n,3,.5)

# generate release lengths and recapture lengths assuming random Linf and measurement error
Linf <- rnorm(n,mu.Linf,sigma.Linf)
sig.tag <- 5
L1 <- Linf*(1-exp(-k*A)) + rnorm(n,0,sig.tag)
L2 <- Linf*(1-exp(-k*(A+dt))) + rnorm(n,0,sig.tag)

# Note that the code allows for a different measurement error for scientist-measured recapture lengths
# versus fisherman-measured recapture lengths, so need to specify whether L2 was measured by
# a scientist (=0) or a fisherman (=1).  Here I'll specify all as scientist (0) since I generated L2
# assuming the same measurement error for all observations.
L2measurer <- rep(0,n)

tagdat <- cbind(L1,L2,rep(0,n),dt,L2measurer)
plot(c(A,A+dt),c(L1,L2),xlab="age",ylab="length")
points(A,L1,col=2)


# b) otolith data:
# ----------------------
set.seed(1)
n.oto <- 200

# generate random uniform ages
age <- runif(n.oto,.5,20)
# generate corresponding lengths with measurement error
Linf <- rnorm(n.oto,mu.Linf,sigma.Linf)
sig.oto <- 5
len <- Linf*(1-exp(-k*(age-a0))) + rnorm(n.oto,0,sig.oto)

otodat <- cbind(age,len)
points(age,len,pch=2,col=3)


# c) length-frequency data:
# -------------------------
set.seed(1)

# Note: the LF component assumes that a model decomposition has already been applied to the
# raw data, such that the dataset input to the growth analysis is the ages and lengths corresponding to
# the modes in the LF data, along with the standard error estimates for each modal length estimate

# generate ages and lengths corresponding to modes in the LF data, assuming we have monthly LF data for
# fish of ages 1-5, and we have this data for 5 years
mode.age <- rep(seq(1,5,1/12),times=5)
n.lf <- length(mode.age)
# recall we are not including a random Linf parameter for the LF data
sig.lf <- 5
mode.len <- mu.Linf*(1-exp(-k*(mode.age-a0))) + rnorm(n.lf,0,sig.lf)

# Specify standard error estimates corresponding to each of the modal length estimates
se.mode <- rlnorm(n.lf,.01,.5)

lfdat <- cbind(mode.age,mode.len,se.mode)
points(mode.age,mode.len,pch=3,col=4)



# Source the R code with the likelihood and a range of growth functions:
# -----------------------------------------------------------------------
source("/Users/stephenscherrer/Google Drive/Weng Lab/Data/Bottomfish/Okamoto's 1990s Mark Recapture Data/src/Laslett Functions/joint_lkhd.r")
source("/Users/stephenscherrer/Google Drive/Weng Lab/Data/Bottomfish/Okamoto's 1990s Mark Recapture Data/src/Laslett Functions/growth_functions.r")



# Fit a VB model WITHOUT seasonal growth:
# -----------------------------------------
growth.ssnl.f<- growthvb.f
npf <- 1  #number of parameters passed to growth.ssnl.f (in this case k)
npA <- 2  #number of parameters in distribution of release ages for tag model

# specify starting parameters, as well as upper and lower bounds
#        mu.L, sig.L,   k,  mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
p0 <- c(  70,   5.0, .10,   1.0,   .10,    5.0,      0,  0,   10.0,   10.0)
lb <- c(   50,   1.0, .05,    .1,   .05,    1.0,      0, -2,    0.0,    0.0)
ub <- c(  100,  15.0, .40,   1.5,   .50,   10.0,      0,  1,   15.0,   15.0)
fit.vb <- nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=tagdat,otodat=otodat,lfdat=lfdatz)
fit.vb


# Fit a VB model *WITH* seasonal growth:
# -----------------------------------------
growth.ssnl.f<- growthvb.ssnl.f
npf <- 3  #number of parameters passed to growth.ssnl.f (in this case k, u, w)
npA <- 2  #number of parameters in distribution of release ages for tag model

# specify starting parameters, as well as upper and lower bounds
#        mu.L, sig.L,   k,  u,   w, mu.A, sig.A, sig.sci, sig.f, a0, sig.oto, sig.lf
p0 <- c(  110,   5.0, .10, .1,   0,  1.0,   .10,    5.0,      0,  0,   10.0,   10.0)
lb <- c(   70,   1.0, .05,  0, -.5,   .1,   .05,    1.0,      0, -2,    0.0,    0.0)
ub <- c(  130,  15.0, .40,  1,  .5,  1.5,   .50,   10.0,      0,  1,   15.0,   15.0)
fit.ssnl.vb <- nlminb(p0,joint.logl.f,lower=lb,upper=ub,npf=npf,npA=npA,tagdat=tagdat,otodat=otodat,lfdat=lfdat)
fit.ssnl.vb




