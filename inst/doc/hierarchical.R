## ----setup, include = FALSE, cache = FALSE, echo = FALSE------------------------------------------
library('ctsem')
#library(knitr)
set.seed(22)
knit_hooks$set(crop = hook_pdfcrop)
opts_chunk$set(warning = FALSE, fig.align = 'center', width.cutoff = 100, fig.show = 'hold', eval = TRUE, echo = TRUE, message = FALSE, comment = NA, tidy = FALSE, out.truncate = 100, size='small', crop=TRUE)
options(width = 100, scipen = 12, digits = 3)


# Tpoints=50
# n.manifest=2
# n.TDpred=1
# n.TIpred=3
# n.latent=2
# n.subjects=20
# gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,n.TDpred=n.TDpred,n.TIpred=n.TIpred,n.manifest=n.manifest,
#   MANIFESTVAR=diag(0.5,2),
#   TIPREDEFFECT=matrix(c(.5,0,0,-.5,0,0),nrow=2),
#   TIPREDVAR=matrix(c(1,-.2,0, 0,1,0, 0,0,.5),nrow=3),
#   TDPREDEFFECT=matrix(c(.1,-.2),nrow=2),
#   TDPREDVAR=matrix(0,nrow=n.TDpred*(Tpoints-1),ncol=n.TDpred*(Tpoints-1)),
#   TDPREDMEANS=matrix(rnorm(n.TDpred*(Tpoints-1),0,1),nrow=n.TDpred*(Tpoints-1)),
#   LAMBDA=diag(1,2), 
#   # DRIFT=matrix(c(-.6+rnorm(1,0,.15),-.2+rnorm(1,0,.1),.12+rnorm(1,0,.1),-.3+rnorm(1,0,.05)),nrow=2),
#   DRIFT=matrix(c(-.3,.2,-.1,-.2),nrow=2),
#   TRAITVAR=t(chol(matrix(c(4,3,3,4),nrow=2))),
#   # T0TRAITEFFECT=diag(3,n.latent),
#   DIFFUSION=matrix(c(.3,.1,0,.2),2),CINT=matrix(c(0,0),nrow=2),T0MEANS=matrix(0,ncol=1,nrow=2),
#   T0VAR=diag(100,2))
# 
# cd<-ctGenerate(gm,n.subjects=n.subjects,burnin=300, dT=1,asymptotes=F,simulTDpredeffect = T)
# model<-ctModel(type='stanct',n.latent=n.latent,n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,LAMBDA=diag(n.latent))
# long<-ctWideToLong(cd,Tpoints,n.manifest=model$n.manifest,manifestNames = model$manifestNames, 
#   n.TDpred=n.TDpred,n.TIpred=n.TIpred,TDpredNames = model$TDpredNames,TIpredNames = model$TIpredNames)
# long<-ctDeintervalise(long)
# long[is.na(long)]<-0
# ctstantestdat <- long


## ----install,eval=FALSE---------------------------------------------------------------------------
#  install.packages("ctsem")

## ----data,echo=FALSE,size='footnotesize'----------------------------------------------------------
temp=ctstantestdat
ctstantestdat[,'time']=ctstantestdat[,'time']+round(rnorm(nrow(ctstantestdat),0,.1),2)
ctstantestdat[c(1:5,87:89),]
ctstantestdat=temp

## ----model----------------------------------------------------------------------------------------
model<-ctModel(type='stanct',
  n.latent=2, latentNames=c('eta1','eta2'),
  n.manifest=2, manifestNames=c('Y1','Y2'),
  n.TDpred=1, TDpredNames='TD1', 
  n.TIpred=3, TIpredNames=c('TI1','TI2','TI3'),
  LAMBDA=diag(2))

## ----modelpars,size='footnotesize'----------------------------------------------------------------
head(model$pars,8)

## ----transform, fig.width=8, fig.height=6---------------------------------------------------------
par(mfrow=c(1,2))
plot(model,rows=11)
print(model$pars$transform[11])
model$pars$transform[11]<- "(exp(param*2) +.0001)*.2"
plot(model,rows=11)

## ----restrictbetween------------------------------------------------------------------------------
model$pars[c(15,18),]$indvarying<-FALSE
model$pars$sdscale[1:28] <- .5

## ----restricttipred-------------------------------------------------------------------------------
model$pars[,c('TI1_effect','TI2_effect','TI3_effect')]<-FALSE
model$pars[c(19,20),c('TI1_effect','TI2_effect','TI3_effect')]<-TRUE

## ----fitting,cache=TRUE---------------------------------------------------------------------------
fit<-ctStanFit(datalong = ctstantestdat, ctstanmodel = model, iter=200, 
  chains=2, plot=FALSE, control=list(max_treedepth = 6))

## ----output,eval=FALSE----------------------------------------------------------------------------
#  summary(fit)

## ----plots1,echo=FALSE,fig.width=8, fig.height=6--------------------------------------------------
ctStanDiscretePars(fit, plot=TRUE)

## ----outputposterior------------------------------------------------------------------------------
ctStanPlotPost(ctstanfitobj = fit, rows=11)

## ----kalmanplot,echo=TRUE,fig.width=10,fig.height=7-----------------------------------------------
ctStanKalman(fit, subjects=2, timerange=c(0,60), timestep=.1, plot=TRUE)

## ----sunspots,cache=TRUE--------------------------------------------------------------------------
#get data
 sunspots<-sunspot.year
 sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
 id <- 1
 time <- 1749:1924
datalong <- cbind(id, time, sunspots)

#setup model
 model <- ctModel(type='stanct', n.latent=2, n.manifest=1, 
  manifestNames='sunspots', 
  latentNames=c('ss_level', 'ss_velocity'),
   LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
   DRIFT=matrix(c(0, 1,   'a21', 'a22'), nrow=2, ncol=2, byrow=TRUE),
   MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
   # MANIFESTVAR=matrix(0, nrow=1, ncol=1),
   CINT=matrix(c(0, 0), nrow=2, ncol=1),
   DIFFUSION=matrix(c(
     0, 0,
     0, "diffusion"), ncol=2, nrow=2, byrow=TRUE))
 
 model$pars$indvarying<-FALSE #Because single subject
 model$pars$transform[14]<- '(param)*5+44 ' #Because not mean centered
 model$pars$transform[4]<-'log(exp(-param*1.5)+1)' #To avoid multi modality 

#fit
fit <- ctStanFit(datalong, model, iter=400, chains=2, control=list(adapt_delta=.9))

#output
summary(fit)$popmeans

## ----sdtransforms---------------------------------------------------------------------------------
#population mean and subject level deviations (pre-transformation)

hypermeans_prior <- rnorm(99999, 0, 1)
hypermeans_post <- -2 #hypothetical sample

indparamsbase_prior <- rnorm(99999, 0, 1)
indparamsbase_post <- rnorm(50, 0, 1) #hypothetical sample

#population standard deviation prior

hypersd_prior <- rnorm(99999, 0, 1)
hypersd_prior <- hypersd_prior[hypersd_prior > 0]

#population standard deviation posterior

hypersd_post <- .4 #hypothetical

#population cholesky correlation matrix
#lower triangle sampled from uniform(-1, 1), 
#upper triangle fixed to 0, 
#diagonal calculated according to hypersd.
hypercorrchol_post <- 1 #because only 1 parameter here...

#population cholesky covariance matrix 
#here based on mean of hypersd_post, for convenience...
#in reality would have multiple samples.
hypercovchol <- diag(hypercorrchol_post,1) %*% 
  diag(hypersd_post,1) %*% diag(hypercorrchol_post,1)

#subject level parameters
#first compute pre transformation parameters
#then transform appropriately (here according to drift auto effect)
indparams <- hypercovchol %*% indparamsbase_post + hypermeans_post 
indparams <- -log(exp(-1.5 * indparams) + 1)

#post transformation population standard deviation
hsd_ourparameter <- abs( #via delta approximation
  (-log(exp(-1.5 * (hypermeans_post + hypersd_post)) + 1) - 
   -log(exp(-1.5 * (hypermeans_post - hypersd_post)) + 1) ) / 2)

## ----priorplots3,fig.height=12,fig.width=10,cache=FALSE,echo=FALSE, crop=FALSE--------------------
param<-rnorm(9999999)
par(mfrow=c(3,2),mar=c(4,3,3,1),mgp=c(2,.5,0),cex=1)
yy<-exp(4*param)
yy<-yy[yy<10]
plot(density(yy,from=-.1 ,bw=.01,n=5000),type='l', lwd=3, xlim=c(-.1,5),ylim=c(0,1),  xaxs='i', yaxs='i',xlab='Value',main='Std. deviation')
grid()
plot(density(yy^2,from=-.1 ,bw=.02,n=5000),type='l', lwd=3, xlim=c(-.5,25), ylim=c(0,.1), xaxs='i', yaxs='i',xlab='Value',main='Variance')
grid()
plot(density(-log(exp(-param*1.5)+1),bw=.02,n=5000),type='l',lwd=3,ylim=c(0,1),xlim=c(-5,0), xaxs='i', yaxs='i', xlab='Value', main='Auto effect')
grid()
plot(density(exp(-log(exp(-param*1.5)+1))),type='l',lwd=3,ylim=c(0,1.6), xaxs='i', yaxs='i', xlab='Value', main='Autoregression | $\\Delta t = 1$')
grid()
plot(density(2/(1+exp(param))-1,n=5000,bw=.03,from=-1,to=1), type='l', lwd=3, xlim=c(-1,1), xaxs='i', ylim=c(0,1), yaxs='i', xlab='Value',main='Partial correlation')
grid()
plot(density((param)),type='l',lwd=2,ylim=c(0,.5),xlim=c(-4,4),xaxs='i',yaxs='i', xlab='Value', main='Other parameters')
grid()

