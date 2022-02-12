
install.packages("locfdr")   ##if locfdr has not been installed

#*****************************************************************************************************
#   This function estimates the sum of heritability explained by all true susceptibility variants in a GWAS, 
#   given the outcome is a continuous trait
#  zall : input the vector of z-statisitics here (note that we recommend the GWAS dataset be pruned first,
#        then the remaining z-statistics can be used as input) 
# totalN : total sample size 
#*****************************************************************************************************

SumVg.cont <- function (zall,totalN) {


D1ss<-function (x, y, xout = x, spar.offset = 0.1384, spl.spar = NULL)
{
    sp <- if (is.null(spl.spar)) {
        sp <- smooth.spline(x, y)
        smooth.spline(x, y, spar = sp$spar + spar.offset)
    }
    else smooth.spline(x, y, spar = spl.spar)
    predict(sp, xout, deriv = 1)$y
}


#***************************************************************************
#                 function to cal. Z from a given Vg			         *
#***************************************************************************

VgtoZ.func <- function(varexp,samp=totalN){
sqrt( varexp/(1-varexp)*(samp-2) )
}


#***************************************************************************
#                 function to cal. Vg from a given z			                 *
#***************************************************************************

ztoVg.func <-function(Z,samp=totalN){
Z^2 /(samp-2 + Z^2)
}


#******************************************************************************
#  Modified Selection bias correction by Efron's empirical Bayes method       *
#******************************************************************************
bias.corr.kernel<- function (zz,xout,bw="nrd0")

{
   density.obj=density(zz,bw=bw)
   fz.func = splinefun(density.obj$x, density.obj$y)
   Psi.z = log(  fz.func(zz) /dnorm(zz)  )
   truez = D1ss(x= zz, y=Psi.z ,xout=xout)
return(truez)
}

#******************************************************************************
#  Selection bias correction by Efron's empirical Bayes method                *
#******************************************************************************
bias.corr <- function (zz, bre = 120, df = 7,xout)

{

      lo <- min(zz)
       up <- max(zz)

    zzz <- pmax(pmin(zz, up), lo)
    breaks <- seq(lo, up, length = bre)
    zh <- hist(zzz, breaks = breaks, plot = F)
    x <- (breaks[-1] + breaks[-length(breaks)])/2  ## x is the mid-point of z-value bins
    yall <- y <- zh$counts

     f <- glm(y ~ ns(x, df = df), poisson)$fit
     
     fz.func = approxfun(x, f, rule = 2, ties = "ordered")
     Psi.z = log(  fz.func(zz) /dnorm(zz)  )
    truez = D1ss(x= zz, y=Psi.z,xout=xout )
return(truez)
}


library(locfdr)


Vgall = ztoVg.func(zall)
zall.corr.ker =  bias.corr.kernel(zall,xout=zall)
Vg.corr.ker = ztoVg.func(zall.corr.ker )
sum.kernel = sum(Vg.corr.ker)

fdrobj2 = locfdr(zall,bre=120,df=7,nulltype=2,plot=0)
fdr2 = fdrobj2$fdr
largefdr.ind = which (fdr2>=0.95) 
truez.given.H1 = zall.corr.ker[-largefdr.ind]/(1-fdr2[-largefdr.ind])   ##note that we don't consider any z-values with fdr>0.95   
Vg.corr = ztoVg.func(truez.given.H1)*(1-fdr2[-largefdr.ind ]) 

sum.kernel.givenH1 = sum(Vg.corr)

return ( list(sum.kernel=sum.kernel, sum.kernel.givenH1=sum.kernel.givenH1) )
}



