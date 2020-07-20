# Negative log-likelihood component for the tag-recapture data
# (using Gauss-Hermite quadrature) and assuming a lognormal distn for A.

# Assumes data is in a matrix with columns: l1, l2, t1, t2, measr
# where l1=release length 
#       l2=recapture length 
#       t1=release time (in julian time since Jan 1 of reference year) 
#       t2=recapture time (in julian time since Jan 1 of reference year)
#       measr=0 if recapture length measured by a scientist
#            =1 if recapture length measured by a fisherman 

# Assumes parameters are: mu.L, sig.L, param.f, mu.logA, sig.logA, sig.s, sig.f
# where param.f is parameters of growth function, mu.logA and sig.logA are
# parameters of the log normal distribution for A, sig.s is the error variance
# for scientists, and sig.f is the error variance for fisherman.

# The log density of the age-at-tagging  must be specified in a function 
# logdensA.f(a,param.A)      
 
# The errors are specified in the function error.ssnl.f

            #------------------------------------------------#


log.lognormA.f<- function(a,param.A)
{   # the log density of the tagging age A 
    # A log-normal 
    alpha.A<- param.A[1]
    beta.A<-  param.A[2]
    
    # a is nf x n matrix
    lna<- a
    tf<- a>0
    
    lna[tf]<- log(a[tf])
    lna[!tf]<- -300

    logpa<- (lna-alpha.A)/beta.A
    logpa<- -0.5*log(2*pi)-log(beta.A)-lna-0.5*logpa^2
    return(logpa)
}

lognormA.f<- function(a,param.A)
{# dlnorm only takes vectors for the first argument
    veca<- as.vector(a)

    alpha.A<- param.A[1]
    beta.A<-  param.A[2]
    
    tf<- a <= 0
    a[tf]<- exp(-30)

    pa<- dlnorm(veca,alpha.A,beta.A)
    pa<- matrix(pa,nrow=nrow(a),ncol=ncol(a),byrow=F)
    return(pa)
}

# The function logdensA.f() is changed if a new distribution for
# A is tried
logdensA.f<- log.lognormA.f
densA.f<-lognormA.f

error.ssnl.f<- function(capstatus,measr,param.sig)
{   # -this gives the additive error structure of the model
    # -capstatus = "release" or "recapture"
    # -measr = 0 if a scientist measured the recapture length
    #       = 1 if a fisherman measured the recapture length
    
    sig.s<- param.sig[1]   # s.d. for a scientiest
    sig.f<- param.sig[2]   	# this is the additional error for fishermen
    if(capstatus=="release")
       {sigma<- sig.s} 
    if(capstatus=="recapture")
       {sigma<- sig.s+sig.f*measr}
    return(sigma)
}


ha.ssnl.f<- function(a,tagdat,param.g,param.sig)
{   # this calculates h(l1,l2|a)
    # it is needed as part of int.f
     
    # a is nf x n matrix
    tf<- a<0
    a[tf]<-0

    h<- exp(logha.ssnl.f(a,tagdat,param.g,param.sig))
    return(h)
}

logha.ssnl.f<- function(a,tagdat,param.g,param.sig)
{   # this calculates log h(l1,l2|a)
    # it is needed as part of logint.f

    # a is nf x n matrix
    tf<- a<0
    a[tf]<-0

    # retrieve the data
    l1<- tagdat[,1]
    l2<- tagdat[,2]
    t1<- tagdat[,3]
    t2<- tagdat[,4]
    measr<- tagdat[,5]
   
    # retrieve the growth curve parameters
    mu.L<-    param.g[1]
    sigma.L<- param.g[2]
    param.f<- param.g[3:length(param.g)]

    ucol<- matrix(1,nrow=nrow(a),ncol=1)
    urow<- matrix(1,nrow=1,ncol=ncol(a))

    a2<- a+(t2-t1)%*%urow
    f1<- growth.ssnl.f(a,param.f,t1-floor(t1))
    f2<- growth.ssnl.f(a2,param.f,t2-floor(t2))
    sigma1<- error.ssnl.f("release",measr,param.sig)
    sigma2<- error.ssnl.f("recapture",measr,param.sig)
    mu1a<- mu.L*f1
    mu2a<- mu.L*f2
    var1a<- (sigma.L*f1)^2 + sigma1^2
    var2a<- (sigma.L*f2)^2 + sigma2^2
    s12<- sqrt(var1a*var2a)
    rhoa<-  (sigma.L^2)*f1*f2/s12
    q1<- l1%*%urow-mu1a
    q2<- l2%*%urow-mu2a
    q12<- q1*q2
    q1<- q1^2/var1a
    q2<- q2^2/var2a
    q12<- q12*rhoa/s12
    q<- q1+q2-2*q12
    logh<- -log(2)-log(pi)-0.5*log(var1a)-0.5*log(var2a)-0.5*log(1-rhoa^2)
    rho<- 2*(1-rhoa^2)
    logh<- logh-q/rho
    return(logh)
}

int.ssnl.f<- function(a,tagdat,param.g,param.A,param.sig)
{   # this calculates h(l1,l2|a) p(a)
    # so that later it can be integrated
    #
    # a is a nf x n matrix, where nf=nrow(tagdat)

    pa<- densA.f(a,param.A)
    integrand<- ha.ssnl.f(a,tagdat,param.g,param.sig)*pa
    return(integrand)
}

logint.ssnl.f<- function(a,tagdat,param.g,param.A,param.sig)
{   # this calculates  -log[h(l1,l2|a) p(a)]
    # it is used in the grid search or golden section search
    # for the mode of h(l1,l2|a) p(a)
    # we search on the log scale because the log is more
    # numerically stable
    #
    # a is a nf x n matrix, where nf=nrow(tagdat)
    # tagdat is a nf x 4 matrix of variables l1,l2,t1,t2

    logpa<- logdensA.f(a,param.A)
    logha<- logha.ssnl.f(a,tagdat,param.g,param.sig)
    integrand<- logha+logpa
    return(integrand)
}

h.ssnl.f<- function(tagdat,param.g,param.A,param.sig)
{
    # This calculates  int_-\infty^\infty h(l1,l2|a) p(a) da
    # using 20-point Gauss-Hermite integration

    gms<- gstats.ssnl.f(tagdat,param.g,param.A,param.sig)
    mng<- gms[,1]
    sdg<- gms[,2]
    pdf.l1.l2<- hermite20.f(int.ssnl.f,tagdat,param.g,param.A,
        param.sig,nrow(tagdat),mng,sdg)
    return(pdf.l1.l2)
}

logl.ssnl.f<-function(param,npf,npA,tagdat)
{   # this calculates the negative log likelihood of the data
    
    npg<- npf+2
    param.g<- param[1:npg]
    param.A<- param[(npg+1):(npg+npA)]
    param.sig<- param[(npg+npA+1):(npg+npA+2)]

    ndata<- nrow(tagdat)
    gms<-   gstats.ssnl.f(tagdat,param.g,param.A,param.sig)
    mng<- gms[,1]
    sdg<- gms[,2]
    pdf.l1.l2<- hermite20.f(int.ssnl.f,tagdat,param.g,param.A,
		param.sig,ndata,mng,sdg)
    neglogl<- -sum(log(pdf.l1.l2))
    if(is.na(neglogl)){print(" neglogl NA ")}
    
    #print(param)
    #print(neglogl)
    return(neglogl)
}


hermite20.f<- function(f,tagdat,param.g,param.A,param.sig,nf,mu,sigma)
 { # this carries out 20-point Gauss-Hermite integration of f(x)
   # i.e. int_-\infty^\infty f(x) dx
   # it assumes that f(x) is roughly Gaussian with mean mu and s.d. sigma
   # here f is a nfx1 vector of functions, and mu, sigma are nfx1 vectors
   #
   # tagdat is a matrix of data which comprise the second argument of f
   # param is a vector of parameters which comprise the third argument of f
   #
   x<- c(-5.3874808900112,-4.6036824495507,-3.9447640401156,-3.3478545673832,
         -2.7888060584281,-2.2549740020893,-1.7385377121166,-1.2340762153953,
         -0.7374737285454,-0.2453407083009, 0.2453407083009, 0.7374737285454,
      1.2340762153953, 1.7385377121166, 2.2549740020893, 2.7888060584281,
      3.3478545673832, 3.9447640401156, 4.6036824495507, 5.3874808900112)

   wexp<- c(0.8985919614532, 0.7043329611769, 0.6222786961914, 0.5752624428525,
            0.5448517423644, 0.5240803509486, 0.5096790271175, 0.4999208713363,
            0.4938433852721, 0.4909215006667, 0.4909215006667, 0.4938433852721,
        0.4999208713363, 0.5096790271175, 0.5240803509486, 0.5448517423644,
        0.5752624428525, 0.6222786961914, 0.7043329611769, 0.8985919614532)

   rt2<- c(1.414213562373)
   urow<- matrix(1,nrow=1,ncol=20)

   xmat<- rt2*(sigma%*%t(x))+mu%*%urow        # xmat is a nf x 20 matrix
   fx<-    f(xmat,tagdat,param.g,param.A,param.sig)
     
   intf<-  fx%*%wexp
   intf<- rt2*sigma*intf
   if(sum(is.na(intf))>0){print(" hermite20 NA ")}
   return(intf)
}

gridsch.f<- function(f,l,u,vmat,param.g,param.A,param.sig)
{   # gridsch.f is a grid search
    #
    # l is vector of true lower bounds
    # u is a vector of true upper bounds

    eps1<- 1e-13

    n<- 100
    urow<- matrix(1,nrow=1,ncol=(n+1))
    range<- u-l
    ind<- (0:n)
    x<- range%*%t(ind)
    x<- l%*%urow+x/n
    fx<- f(x,vmat,param.g,param.A,param.sig)
    fmax<- apply(fx,1,max)
    xmax<- fmax%*%urow
    tf<- fx >= xmax - eps1

    # if there are two minima in a row by some extraordinary chance,
    # xmax found in the following line of code will keep both
    # this is undesirable
 
    # xmax<- as.vector(t(x)[t(tf)])

    # the following piece of code avoids that problem
    rmax<- row(x)[tf]
    cmax<- col(x)[tf]
    tf<- (!duplicated(rmax))
    rmax<- rmax[tf]
    cmax<- cmax[tf]
    ord<- order(rmax)
    rmax<- rmax[ord]
    cmax<- cmax[ord]
    rcmax<- cbind(rmax,cmax)
    xmax<- as.vector(x[rcmax])
    if(min(cmax)==1) {print("maximum on lower bound")
                          print("rows")
                          print(rmax[cmax==1])
                          #browser()
                          cmax[cmax==1]<-2
    }
    if(max(cmax)==ncol(x)) {print("maximum on upper bound")
                          print("rows")
                          print(rmax[cmax==ncol(x)])
                          #browser()
                          cmax[cmax==ncol(x)]<-ncol(x)-1
    }
    rcmax<- cbind(rmax,cmax-1)
    a1<- as.vector(x[rcmax])
    f1<- as.vector(fx[rcmax])
    rcmax<- cbind(rmax,cmax+1)
    a3<- as.vector(x[rcmax])
    f3<- as.vector(fx[rcmax])
    return(cbind(a1,f1,xmax,fmax,a3,f3))
}

quadmax.f<- function(ahdat)
{   # we are attempting to improve the estimate of the location 
    # of the maximum of a nearly quadratic function h(a)
    # ahdat = cbind(a1,h1,a2,h2,a3,h3)
    # a1 < a2 < a3
    # h2>h1 and h2>h3

    a1<- ahdat[,1]
    h1<- ahdat[,2]
    a2<- ahdat[,3]
    h2<- ahdat[,4]
    a3<- ahdat[,5]
    h3<- ahdat[,6]

    dett<-  (a1-a2)*(a3-a2)*(a3-a1)
    alpha<- ((a3-a2)^2)*(h1-h2)-((a1-a2)^2)*(h3-h2)
    beta<-  -(a3-a2)*(h1-h2)+(a1-a2)*(h3-h2)
    alpha<- alpha/dett
    beta<-  beta/dett
    mng<- a2-alpha/(2*beta)
    sdg<- sqrt(-1/(2*beta))
    return(cbind(mng,sdg))
}


gstats.ssnl.f<- function(tagdat,param.g,param.A,param.sig)
{   qt1<- 0.001
    qt2<- 15
    
    l<- matrix(rep(qt1,nrow(tagdat)),ncol=1)
    u<- matrix(rep(qt2,nrow(tagdat)),ncol=1)
   
    aa<- gridsch.f(logint.ssnl.f,l,u,tagdat,param.g,param.A,param.sig)
    brk<-0
    for(i in (1:10))
       {ms<- quadmax.f(aa)

	# we need to make sure that the mng values have stabilised
	# if it has we break out
	if(i>3)
	   {stdiff<- abs((ms[,1]-mngold)/ms[,2])
	    msd<- max(stdiff)
	    if(msd<0.2)
	    {   brk<- 1
	        break
	       }
           }
	
	# if not, we refine the estimates
        a2<- pmax(ms[,1],0.9*aa[,1]+0.1*aa[,5])
        a2<- pmin(a2,0.1*aa[,1]+0.9*aa[,5])
        a1<- pmin(a2,aa[,3])
        a2<- a2+aa[,3]-a1
        a1<- pmin(a1,0.1*aa[,1]+0.9*a2)
        f12<- logint.ssnl.f(cbind(a1,a2),tagdat,param.g,param.A,param.sig)
        tf<- f12[,1]<f12[,2]
        atemp<- cbind(a1,f12[,1],a2,f12[,2],aa[,5:6])
        aa<- cbind(aa[,1:2],atemp[,1:4])
        aa[tf,]<-atemp[tf,]
	mngold<- ms[,1]
       }
    if(brk==0) {assign("nits.gstats",1,where=1)
		print("max # iterations reached")}
    ms<- quadmax.f(aa)
    return(ms)
}


a.l1.l2.ssnl.f<- function(a,tagdat,param.g,param.A,param.sig)
{   # This calculates  a h(l1,l2|a) p(a)

    apa<- a*densA.f(a,param.A)
    integrand<- ha.ssnl.f(a,tagdat,param.g,param.sig)*apa
    return(integrand)
}

ahat.ssnl.f<- function(tagdat,param.g,param.A,param.sig)
{   # this estimates E[A|l1,l2] from the data
    # for given parameter values

    gms<-   gstats.ssnl.f(tagdat,param.g,param.A,param.sig)
    mng<- gms[,1]
    sdg<- gms[,2]
    mean.A<- hermite20.f(a.l1.l2.ssnl.f,tagdat,param.g,param.A,param.sig,
		nrow(tagdat),mng,sdg)
    pdf.l1.l2<- h.ssnl.f(tagdat,param.g,param.A,param.sig)
    mean.A<- mean.A/pdf.l1.l2
    return(mean.A)
}


linfhat.a.ssnl.f<- function(a,tagdat,param.g,param.A,param.sig)
{	# this calculates E[linf|l1,l2,a]
	
	l1<-	tagdat[,1]
	l2<- 	tagdat[,2]
	t1<-	tagdat[,3]
	t2<- 	tagdat[,4]
	measr<- tagdat[,5]

	mu.L<-		param.g[1]
	sig.L<-	param.g[2]
	param.f<- param.g[3:length(param.g)]
	
	ucol<- matrix(1,nrow=nrow(a),ncol=1)
  	urow<- matrix(1,nrow=1,ncol=ncol(a))
  	a2<- a+(t2-t1)%*%urow
  	f1<- growth.ssnl.f(a,param.f,t1-floor(t1))
  	f2<- growth.ssnl.f(a2,param.f,t2-floor(t2))
	sig1<- 	error.ssnl.f("release",measr,param.sig)
   	sig2<- 	error.ssnl.f("recapture",measr,param.sig)

	denom<- sig1^2*sig2^2 + sig.L^2*(sig1^2*f2^2 + sig2^2*f1^2)	
	numer<- sig2^2*f1*(l1 - mu.L*f1) + sig1^2*f2*(l2 - mu.L*f2)
	linfhat.a<- mu.L + sig.L^2*numer/denom

	return(linfhat.a)
}

linf.int.ssnl.f<- function(a,tagdat,param.g,param.A,param.sig)
{	# this calculates E[linf|l1,l2,a]*h(l1,l2|a)*p(a)

	integrand<- linfhat.a.ssnl.f(a,tagdat,param.g,param.A,param.sig)*
		ha.ssnl.f(a,tagdat,param.g,param.sig)*densA.f(a,param.A)
}

linfhat.ssnl.f<-function(tagdat,param.g,param.A,param.sig)
{ 	# this estimates E[Linf|l1,l2] from the data
  	# for given parameter values

	gms<-   gstats.ssnl.f(tagdat,param.g,param.A,param.sig)
   	mng<- gms[,1]
   	sdg<- gms[,2]
	mean.linf<- hermite20.f(linf.int.ssnl.f,tagdat,param.g,param.A,
		param.sig,nrow(tagdat),mng,sdg) 
	pdf.l1.l2<- h.ssnl.f(tagdat,param.g,param.A,param.sig)
	mean.linf<- mean.linf/pdf.l1.l2
   return(mean.linf)
}


