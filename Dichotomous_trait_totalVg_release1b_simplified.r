zvec =      ##MUST enter: the vector of z-values from which total Vg will be estimated
                  ## we advise LD pruning to be performed beforehand

fdrtrunc = 0.95    ##threshold to truncate local fdr values
repl = 10 ##number of simulations to be performed [ie no. of times to assign random signs to z-values]
K =0.001  ##overall disease probability
caseNo = 3230   ##no. of cases
ctrlNo = 4829   ##no. of controls

library(locfdr)
library(sfsmisc)
library(numDeriv)

muaa=0 ##fixed
preval = K 
#set.seed(333)

sum.kernel = numeric(repl)
sum.kernel.givenH1 = numeric(repl)

D1ss<-function (x, y, xout = x, spar.offset = 0.1384, spl.spar = NULL)
{
    sp <- if (is.null(spl.spar)) {
        sp <- smooth.spline(x, y)
        smooth.spline(x, y, spar = sp$spar + spar.offset)
    }
    else smooth.spline(x, y, spar = spl.spar)
    predict(sp, xout, deriv = 1)$y
}

func.RR <- function(RR1,Vg,PA=0.5){
RR2=RR1^2
Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2

faa = K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa 
T = qnorm(1-faa) #muaa is set to 0 and residual var set to 1
muAa = T-qnorm(1-fAa)
muAA = T-qnorm(1-fAA)

mean.all= PAa*muAa+ PAA*muAA
expVg= Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
expVg2 = expVg/(1+expVg)
return( (expVg2-Vg)^2 ) 
}

####function to cal. the resulting RR from a given Vg##########
resRR.func <-function(varexp){
optimize(func.RR,c(1,1000),Vg=varexp)$minimum}


##function to cal. power ##
power.func <- function(RR1,RR2=RR1^2, PA=0.5, K, case.size=caseNo,ctrl.size=ctrlNo,alpha=5e-05){
Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2
muaa=0 #fixed
faa = K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa 

Paa.case = faa*Paa/K
PAa.case = fAa*PAa/K
PAA.case = fAA*PAA/K

Paa.ctrl = (1-faa)*Paa/(1-K)
PAa.ctrl = (1-fAa)*PAa/(1-K)
PAA.ctrl = (1-fAA)*PAA/(1-K)

PA.case = PAa.case/2 + PAA.case
Pa.case = 1- PA.case
PA.ctrl = PAa.ctrl/2 + PAA.ctrl 
Pa.ctrl = 1- PA.ctrl

ORall = PA.case*Pa.ctrl/Pa.case/PA.ctrl
VarlnOR = 1/(2*case.size)*(1/PA.case + 1/Pa.case) +  1/(2*ctrl.size)*(1/PA.ctrl +1/Pa.ctrl)
Z = log(ORall)/sqrt(VarlnOR)


critR = qnorm(1-alpha/2)
critL = qnorm(alpha/2)

power = 1-pnorm(critR,mean=Z,sd=1)+pnorm(critL,mean=Z,sd=1)
res.power.func <- list()
res.power.func$Z = Z
res.power.func$power = power
return(res.power.func) 
}

muaa=0 #fixed

##function to be optimized##
func.RR <- function(RR1,Vg,PA=0.5){
RR2=RR1^2
Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2

faa = K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa 
T = qnorm(1-faa) #muaa is set to 0 and residual var set to 1
muAa = T-qnorm(1-fAa)
muAA = T-qnorm(1-fAA)
##overall mean
mean.all= PAa*muAa+ PAA*muAA
expVg= Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
expVg2 = expVg/(1+expVg)
return( (expVg2-Vg)^2 ) 
}

#***************************************************************************
#                 function to cal. Z from a given Vg			         *
#***************************************************************************

VgtoZ.func <- function(varexp,PA=0.5, K=0.001, case.size=caseNo,ctrl.size=ctrlNo){

RR1 = optimize(func.RR,c(1,100000),Vg=varexp)$minimum

RR2=RR1^2
Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2
muaa=0 #fixed
faa = K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa 

Paa.case = faa*Paa/K
PAa.case = fAa*PAa/K
PAA.case = fAA*PAA/K

Paa.ctrl = (1-faa)*Paa/(1-K)
PAa.ctrl = (1-fAa)*PAa/(1-K)
PAA.ctrl = (1-fAA)*PAA/(1-K)

PA.case = PAa.case/2 + PAA.case
Pa.case = 1- PA.case
PA.ctrl = PAa.ctrl/2 + PAA.ctrl 
Pa.ctrl = 1- PA.ctrl

Acase = 2*case.size*PA.case
acase = 2*case.size*Pa.case
Actrl = 2*ctrl.size*PA.ctrl
actrl = 2*ctrl.size*Pa.ctrl
chisqmat = matrix(c(Acase,acase,Actrl,actrl),nrow=2, byrow=F)
Zsq = chisq.test(chisqmat)$statistic
Z = sqrt(Zsq)
return(Z) 
}


##function to be optimized
power.optim <- function(RR1,obsZ,RR2=RR1^2,case.size=caseNo,ctrl.size=ctrlNo,PA=0.5, K=preval){
Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2
muaa=0 #fixed
faa = K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa

Paa.case = faa*Paa/K
PAa.case = fAa*PAa/K
PAA.case = fAA*PAA/K

Paa.ctrl = (1-faa)*Paa/(1-K)
PAa.ctrl = (1-fAa)*PAa/(1-K)
PAA.ctrl = (1-fAA)*PAA/(1-K)

PA.case = PAa.case/2 + PAA.case
Pa.case = 1- PA.case
PA.ctrl = PAa.ctrl/2 + PAA.ctrl
Pa.ctrl = 1- PA.ctrl

ORall = PA.case*Pa.ctrl/Pa.case/PA.ctrl
							
VarlnOR = 1/(2*case.size)*(1/PA.case + 1/Pa.case) +  1/(2*ctrl.size)*(1/PA.ctrl +1/Pa.ctrl)
Z = log(ORall)/sqrt(VarlnOR)

return( (Z-obsZ)^2 )
}

##########################################################################################
ztoVg.func <-function(obsZ,PA=0.5,K=preval){

RR1=nlminb(1.13,power.optim,obsZ=obsZ,case.size=caseNo,ctrl.size=ctrlNo,PA=0.5,K=preval)$par
RR2=RR1^2

Paa = (1-PA)^2
PAa = 2*PA*(1-PA)
PAA = PA^2
muaa=0 #fixed

faa= K/(Paa + PAa*RR1 + PAA*RR2)
fAa= RR1*faa
fAA= RR2*faa

T = qnorm(1-faa) #muaa is set to 0 and residual var set to 1
muAa = T-qnorm(1-fAa)
muAA = T-qnorm(1-fAA)

##overall mean
mean.all= PAa*muAa+ PAA*muAA
Vg= Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
return( Vg/(1+Vg) )
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

NoOfMarkers=length(zvec)

t1=proc.time()
#******************************************************************
#              start of simulation loop                           *
#******************************************************************
for (j in 1:repl) {

randno = runif( NoOfMarkers , -1,1) 
zall=sign(randno)* zvec

zall.corr.ker =  bias.corr.kernel(zall,xout=zall)
Vg.corr.ker = mapply( ztoVg.func, obsZ= zall.corr.ker  )
sum.kernel[j] = sum(Vg.corr.ker)

fdr = locfdr(zall,bre=500,df=10,nulltype=2,plot=0)$fdr
largefdr.ind = which (fdr>=fdrtrunc) 
truez.given.H1 = zall.corr.ker[-largefdr.ind]/(1-fdr[-largefdr.ind])   ##note that we don't consider any z-values with fdr>0.95   
Vg.corr = mapply( ztoVg.func, obsZ= truez.given.H1 )*(1-fdr[-largefdr.ind]) 
sum.kernel.givenH1[j] = sum(Vg.corr)
}  #***********************************end of simulation loop ****************************
proc.time()-t1

##The results 
res = cbind(sum.kernel,sum.kernel.givenH1)
print( res ) 

## Mean from averaging multiple runs (each time with random signs assigned to z-values)
print( apply(res,2,mean) ) 


