
# In the integrated analysis, all of the data components share the same mean Linf parameter ("mu.Linf").
# The tagging and otolith components include a random Linf parameter, so also share the same "sig.Linf";
# whereas the length-frequency component does not include a random Linf (thus no "sig.Linf") since
# the data are for young fish that are not near their asymptotic length.

# There is no information about the a0 parameter from the tagging data, so this parameter is estimated
# from the otolith and length-frequency data sets

# The variance parameters are assumed to be unique for each data set

# The length-frequency (LF) component assumes that a model decomposition has already been applied to the
# raw data, such that the dataset input to the growth analysis is the lengths and ages corresponding to
# the modes in the LF data.

# - param.f contains the parameters of the growth function, f (note L=Linf*f, so param.f does not include Linf)
#   (e.g. param.f={k} for the standard VB model; param.f={k,u,w} for VB with seasonal growth)
# - npf = number of parameters for the growth function, f; i.e., length of param.f
#   (e.g. npf=1 for standard VB; npf=3 for VB with seasonal growth)
# - param.A contains the parameters for the distribution of A (e.g. log normal in our example)
# - npA = number of parameters for distribution of A, i.e., length of param.A (e.g. npA=2 for log normal)
# - param contains all of the parameters being estimated, including the variance parameters; note that
#   for the tagging component, we allow measurement error in length to be different for scientists vs fishermen
#   since many of our recaptures were measured by fishermen, hence the 2 variance parameters sig.s and sig.f
#   (e.g., param = {mu.Linf, sig.Linf, param.f, param.A, sig.s, sig.f, a0, sig.oto, sig.lf}


# 1) Otolith negative log-likelihood component:
# - assumes data is a matrix with columns: age, length

logl.oto.f <- function(param, npf, npA, otodat)
{
  #print(param)
  # Get data:
	age <- otodat[, 1]
	len <- otodat[, 2]
  n <- nrow(otodat)

	# Get param:
	mu.L <- param[1]
	sig.L <- param[2]
	param.f <- param[3:(npf+2)]
	a0 <- param[2+npf+npA+3]
	sig.oto <- param[2+npf+npA+4]
	
	grwth <- growth.ssnl.f(age - a0, param.f, age)
	exp.len <- mu.L * grwth
	var.len <- (grwth * sig.L)^2 + sig.oto^2
	neglogl <- 1/2 * sum(log(2 * pi * var.len) + ((len - exp.len)^2)/var.len)
	
	#print(neglogl)
	return(neglogl)
}


#2) Length-frequency (LF) negative log-likelihood component:
# - assumes data is in a matrix with columns: age, mean.mode, se.mode

logl.lf.f <- function(param, npf, npA, lfdat)
{
	# Get data:
	age <- lfdat[, 1]
	mu.mode <- lfdat[, 2]
	se.mode <- lfdat[, 3]
	n <- nrow(lfdat)

	# Get parameters:
	mu.L <- param[1]
	param.f <- param[3:(npf+2)]
	a0 <- param[2+npf+npA+3]
	#a0 = 0 # Note: Uncomment this to remove the effect of fitting a0.
	sig.lf <- param[2+npf+npA+5]

	grwth <- growth.ssnl.f(age - a0, param.f, age)
	exp.mode <- mu.L * grwth
	var.mode <- sig.lf^2 + se.mode^2
	neglogl <- 1/2 * sum(log(2 * pi * var.mode) + ((mu.mode - exp.mode)^2)/var.mode)
	#print(param)
	#print(neglogl)
	return(neglogl)
}

#3) Tagging negative log-likelihood component:
# - code and description found in file 'tag lkhd.r'
#source("/Users/stephenscherrer/Google Drive/Weng Lab/Data/Bottomfish/Okamoto's 1990s Mark Recapture Data/src/Laslett Functions/tag_lkhd.r")


# JOINT LIKELIHOOD:
joint.logl.f <- function(param,npf,npA,tagdat,otodat,lfdat,wt.oto=1,wt.tag=1,wt.lf=1)
{
  neglogl.tag<- 0
  neglogl.oto<- 0
  neglogl.lf<- 0
  if(wt.tag>0)  neglogl.tag <- logl.ssnl.f(param,npf,npA,tagdat)
  if(wt.oto>0)  neglogl.oto <- logl.oto.f(param,npf,npA,otodat)
  if(wt.lf>0) neglogl.lf <- logl.lf.f(param,npf,npA,lfdat)
  neglogl <- wt.tag*neglogl.tag + wt.oto*neglogl.oto + wt.lf*neglogl.lf
  #print(param)
  #print(c(neglogl.tag, neglogl.oto, neglogl.lf, neglogl))
  return(neglogl)
}

# where...

# logl.ssnl.f is for tagging data and is found in "tag lkhd.r" (along with all
# necessary associated functions)




