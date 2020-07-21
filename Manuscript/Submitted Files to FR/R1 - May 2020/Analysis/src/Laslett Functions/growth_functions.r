
# Possible growth functions:
# ---------------------------

# L(t) = Linf * f(t-t0)

# where f(t-t0) can be any of the following functions:


growthvb.f<- function(a,param.f,ssnl=0)
{   # this is the standard VB growth curve
    tf<- a<0
    a[tf]<- 0
    k<- param.f[1]
    f<- 1-exp(-k*a)
    return(f)
}


growthvblogk.f<- function(a,param.f,ssnl=0)
{   # this is the VB log k growth curve
    tf<- a<0
    a[tf]<- 0
    k1<-    param.f[1]
    k2<-    param.f[2]
    alpha<- param.f[3]
    beta<-  param.f[4]
    theta<- -(k2-k1)/beta
    d1<- 1+exp(-beta*(a-alpha))
    d2<- 1+exp(beta*alpha)
    g<-  (d1/d2)^theta
    g<-  1-exp(-k2*a)*g
    return(g)
}

growthrich.f <- function(a, param.f,ssnl=0)
{	# this is the Richard's growth curve
	#tf<- a<0
	#a[tf]<- 0
	k <- param.f[1]
	alpha <- param.f[2]
	beta <- param.f[3]
	k <- param.f[1]
	alpha <- param.f[2]
	beta <- param.f[3]
	g <- (1 + beta * exp( - k * (a - alpha)))^(-1/beta) - (1 + beta * exp(k * alpha))^
		(-1/beta)
	return(g)
}


# Seasonal growth models
#--------------------------------------------------------

growthvb.ssnl.f<- function(a,param.f,julian)
{   # this is the VB growth curve with a sinusoidal seasonal component
    # "julian" is in decimal yrs (ie. fraction of yr since Jan 1)
    
    tf<- a<0
    a[tf]<- 0
    
    k<- param.f[1]
    u<- param.f[2]
    w<- param.f[3]
    1-exp(-k*(a+u/(2*pi)*sin(2*pi*(julian-w))))
}


growthvblogk.ssnl.f<- function(a,param.f,julian)
{   # this is the VB log k growth curve
    # "julian" is in decimal yrs (ie fraction of yr since Jan 1)
    
    tf<- a<0
    a[tf]<- 0
    
    k1<-    param.f[1]
    k2<-    param.f[2]
    alpha<- param.f[3]
    beta<-  param.f[4]
    u<-     param.f[5]
    w<-     param.f[6]

    theta<- -(k2-k1)/beta
    d1<- 1+exp(-beta*(a+u/(2*pi)*sin(2*pi*(julian-w))-alpha))
    d2<- 1+exp(beta*alpha)
    g<-  (d1/d2)^theta
    g<-  1-exp(-k2*(a+u/(2*pi)*sin(2*pi*(julian-w))))*g
    return(g)
}









