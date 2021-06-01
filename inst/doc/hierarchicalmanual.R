## ----setup, include = FALSE, cache = FALSE, echo = FALSE------------------------------------------
library('ctsem')
library(knitr)
set.seed(22)
knit_hooks$set(crop = hook_pdfcrop)
opts_chunk$set(warning = FALSE, fig.align = 'center', width.cutoff = 100, fig.show = 'hold', eval = TRUE, echo = TRUE, message = FALSE, comment = NA, tidy = FALSE, autodep=TRUE, out.truncate = 100, size='small', crop=TRUE, fig.pos="htbp",pardef=TRUE,cache=TRUE)
knit_hooks$set(pardef = function(before, options, envir) {
    if (before) par(mfrow=c(1,1),mgp=c(1.5,.6,0),mar=c(3,2,2,1)+.2, cex=.7)  
})
options(width = 100, scipen = 12, digits = 3)

set.seed(1)

library(ggplot2)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}





## ----install,eval=FALSE---------------------------------------------------------------------------
#  install.packages("ctsem")
#  library("ctsem")

## ----data,echo=FALSE,size='footnotesize'----------------------------------------------------------
ctstantestdat[c(3:5,17:22),]

## ----scaling,echo=TRUE----------------------------------------------------------------------------
ctstantestdat[,c('Y1','Y2','TI1','TI2','TI3')] <- 
 scale(ctstantestdat[,c('Y1','Y2','TI1','TI2','TI3')])

## ----model----------------------------------------------------------------------------------------
model<-ctModel(type='stanct',
  latentNames=c('eta1','eta2'),
  manifestNames=c('Y1','Y2'),
  CINT=matrix(c('cint1','cint2 |||| TI1'),nrow=2,ncol=1),
  MANIFESTMEANS=matrix(c(0,0),nrow=2,ncol=1),
  TDpredNames='TD1', 
  TIpredNames = 'TI1',tipredDefault = FALSE,
  LAMBDA=diag(2)) 

## ----model2---------------------------------------------------------------------------------------
model<-ctModel(type='stanct',
  latentNames=c('eta1','eta2'),
  manifestNames=c('Y1','Y2'),
  CINT=c('cint1','cint2 |||| TI1'),
  MANIFESTMEANS=0,
  TDpredNames='TD1', 
  TIpredNames = 'TI1',tipredDefault = FALSE,
  LAMBDA=c(1,0,0,1))

## ----texout,results='asis',echo=FALSE-------------------------------------------------------------
cat(paste0('\\begin{sidewaysfigure}'))
cat(ctModelLatex(model,textsize = 'tiny',compile = FALSE,equationonly = TRUE))
cat('\\caption{Output from ctModelLatex function }
\\label{fig:ctmodellatexout}
\\end{sidewaysfigure}') 

## ----modelpars,size='footnotesize'----------------------------------------------------------------
head(model$pars,4)

## ----fitting,include=FALSE,cache=TRUE-------------------------------------------------------------
if(w32chk()){
suppressWarnings(fit<-ctStanFit(datalong = ctstantestdat, 
  ctstanmodel = model,optimize=TRUE,cores=1,
  # savescores=TRUE,savesubjectmatrices = TRUE, 
  nopriors=FALSE))
}

## ----fittingshow,include=TRUE,eval=FALSE----------------------------------------------------------
#  fit<-ctStanFit(datalong = ctstantestdat, ctstanmodel = model, nopriors=FALSE)

## ----output,eval=FALSE----------------------------------------------------------------------------
#  summary(fit, timeinterval = 1)

## ----ctStanContinuousPars,eval=FALSE--------------------------------------------------------------
#  ctStanContinuousPars(fit,calcfunc = quantile, calcfuncargs = list(probs=.975))

## ----plots1,echo=2,fig.width=5, fig.height=3,fig.cap='Discrete-time cross-effect dynamics of the estimated system for a range of time intervals, with 95\\% credible intervals.'----
if(w32chk()){try({
ctStanDiscretePars(fit, plot=TRUE, indices = 'CR')
})}

## ----outputposterior, echo=2,fig.width=6, fig.height=3, fig.cap='Prior and posterior densities relevant to the first variables manifest intercept.'----
if(w32chk()){
ctStanPlotPost(obj = fit, rows=3) 
}

## ----kalmanplot, echo=c(-1,-4),fig.width=8, fig.height=4, fig.cap='Predictions for three subjects over two processes. Uncertainty shown is a 95\\% credible interval comprising both process and measurement error.'----
if(w32chk()){
ctKalman(fit, subjects=c(2,4,5),  kalmanvec=c('y', 'yprior'), 
  plot=TRUE, timestep=.01)
}

## ----tipredeffects,echo=c(2,3), fig.width=5, fig.height=4, fig.cap='Expectations for individuals parameter values change depending on their score on time independent predictors.'----
if(w32chk()){
ctStanTIpredeffects(fit, plot = TRUE, whichpars=c('dtDRIFT[2,1]','CINT'), 
 timeinterval = .5, whichTIpreds = 1, includeMeanUncertainty = TRUE) 
}

## ----echo=FALSE,transform,fig.align='center', fig.height=4,fig.cap="Prior distribution density plots."----
#Manual change after model specification
if(w32chk()){
model$pars$indvarying[7] <- TRUE
p <- plot(model, rows=7,rawpopsd=1, plot=FALSE)[[1]] + 
  ggplot2::theme(legend.position = "none")+ggplot2::ggtitle('Default prior')

model$pars$transform[7]<- '2 * param -1'

p2 <- plot(model, rows=7,rawpopsd=1,plot=FALSE)[[1]]

p3 <- g_legend(p2)
p2 <- p2 + ggplot2::theme(legend.position ='none')+ggplot2::ggtitle('Modified prior')
gridExtra::grid.arrange(p,p2,p3,layout_matrix=matrix(c(3,1,1,1,1,1,1,1,3,2,2,2,2,2,2,2),ncol=2))
}

## ----echo=TRUE,eval=FALSE,transformne-------------------------------------------------------------
#  #Manual change after model specification
#  model$pars$indvarying[7] <- TRUE
#  p <- plot(model, rows=7,rawpopsd=1, plot=FALSE)
#  print(p[[1]] + ggplot2::theme(legend.position = "none"))
#  model$pars$transform[7]<- '2 * param -1'
#  
#  #Change during model specification
#  model<-ctModel(type='stanct',
#  	DRIFT = matrix(c('drift_eta1_eta1 | 2 * param -1 | TRUE',
#  	  'dr21','dr12','dr22'),2,2),
#  	latentNames=c('eta1','eta2'),
#  	manifestNames=c('Y1','Y2'),
#  	TDpredNames='TD1',
#  	TIpredNames=c('TI1'),
#  	LAMBDA=diag(2))
#  
#  plot(model, rows=7,rawpopsd=1)

## ----transform2,  echo=FALSE,fig.height=4,fig.align='center',fig.cap="Prior distribution density plots of auto-effects, with default (left) and adjusted (right) scale parameter for population standard deviation. "----
if(w32chk()){
p <-plot(model, rows=7, plot=FALSE,nsamples = 100000)[[1]] + ggplot2::theme(legend.position ='none')
model$pars$sdscale<- 3
p2 <- plot(model, rows=7,plot=FALSE,nsamples=100000)[[1]]
p3 <- g_legend(p2)
p2 <- p2 + ggplot2::theme(legend.position ='none')
g=list(p,p2) 
lim <- sapply(g,function(p){
  g=ggplot_build(p)$layout$panel_params[[1]]
  return(c(x=g$x.range,y=g$y.range))
})
g <- lapply(g,function(x) x + coord_cartesian(xlim=c(-5,2.5), 
  ylim=c(mean(lim[3,]),mean(lim[4,]) )))
g[[1]] <- g[[1]] +ggtitle('Default sdscale (1.0)')
g[[2]] <- g[[2]] +ggtitle('Modified sdscale (3.0)')
gridExtra::grid.arrange(g[[1]],g[[2]],p3,layout_matrix=matrix(c(3,1,1,1,1,1,1,1,3,2,2,2,2,2,2,2),ncol=2))
}

## ----transform2ne,  eval=FALSE--------------------------------------------------------------------
#  plot(model, rows=7)
#  model$pars$sdscale<- .1
#  plot(model, rows=7)

## ----restrictbetween------------------------------------------------------------------------------
model$pars$indvarying[ !(model$pars$matrix %in% c('DRIFT','MANIFESTMEANS'))] <- FALSE
model$pars$indvarying[ (model$pars$matrix %in% c('DRIFT','MANIFESTMEANS'))] <- TRUE

## ----restricttipred-------------------------------------------------------------------------------
model$pars[,c('TI1_effect','TI2_effect','TI3_effect')] <- FALSE
model$pars[model$pars$matrix == 'DRIFT', 
  c('TI1_effect','TI2_effect','TI3_effect')] <- TRUE

## ----sunspots,include=FALSE, cache=TRUE-----------------------------------------------------------
if(w32chk()){
#get data
 sunspots<-sunspot.year
 sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
 id <- 1
 time <- 1749:1924
datalong <- cbind(id, time, sunspots)

#setup model
 ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1, 
  manifestNames='sunspots', 
  latentNames=c('ss_level', 'ss_velocity'),
   LAMBDA=matrix(c( 1, 'ma1| log(1+(exp(param)))' ), nrow=1, ncol=2),
   DRIFT=matrix(c(0, 'a21 | -log(1+exp(param))', 1, 'a22'), nrow=2, ncol=2),
   MANIFESTMEANS=matrix(c('m1|param * 10 + 44'), nrow=1, ncol=1),
   MANIFESTVAR=diag(0,1), #As per original spec
   CINT=matrix(c(0, 0), nrow=2, ncol=1),
   DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))

#fit
ssfit <- ctStanFit(datalong, ssmodel, 
  iter=300, chains=2,optimize=FALSE, nopriors=FALSE)

#output
summary(ssfit)$popmeans
}

## ----sunspotsshow,include=TRUE,eval=FALSE---------------------------------------------------------
#  #get data
#  sunspots<-sunspot.year
#  sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
#  id <- 1
#  time <- 1749:1924
#  datalong <- cbind(id, time, sunspots)
#  
#  #setup model
#  ssmodel <- ctModel(type='stanct',
#    manifestNames='sunspots',
#    latentNames=c('ss_level', 'ss_velocity'),
#    LAMBDA=c( 1, 'ma1| log(1+(exp(param)))'),
#    DRIFT=c(0, 1,
#      'a21 | -log(1+exp(param))','a22'),
#    MANIFESTMEANS=c('m1|param * 10 + 44'),
#    MANIFESTVAR=diag(0,1), #As per original spec
#    CINT=0,
#    DIFFUSION=c(0, 0,
#      0, "diffusion"))
#  #fit
#  ssfit <- ctStanFit(datalong, ssmodel,
#    iter=300, chains=2, optimize=FALSE, nopriors=FALSE)
#  
#  #output
#  summary(ssfit)$popmeans
#  

## ----nl,eval=FALSE,results='hide',messages=FALSE,cache=TRUE---------------------------------------
#  sunspots<-sunspot.year
#  sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
#  id <- 1
#  time <- 1749:1924
#  datalong <- data.frame(id, time, sunspots)
#  
#  m <- ctModel(type='stanct',
#    manifestNames='sunspots',
#    latentNames=c('ss_level', 'ss_velocity'),
#    LAMBDA=c( 1, 'ma1|log(1+exp(param))'),
#    DRIFT=c(0, 1,
#      '-log1p_exp(freqintercept + freqbylevel *  ss_level)','a22'),
#    MANIFESTMEANS=c('m1|param * 10 + 44'),
#    MANIFESTVAR=diag(0,1), #As per original spec
#    CINT=0,
#    DIFFUSION=c(0, 0,
#      0, "diffusion"),
#    PARS=c('freqintercept', 'freqbylevel'))
#  
#  ssfitnl <- ctStanFit(datalong, m)

## ----transforms,  fig.width=8, fig.height=6, fig.cap='Depiction of the prior distributions and sampling process through which individual specific parameters are determined. Note that the population sd is strictly positive, however the density plots involve some smoothing."'----
#set plotting parameters
par(mfrow=c(2,2), lwd=3, yaxs='i', mgp=c(1.8,.5,0), 
  mar=c(3,3,3,1)+.1)
bw=.03 

n <- 999999 #number of samples to draw to from prior for plotting purposes
nsubjects <- 4 #number of subjects

#parameter specific transform
tform <- function(x) -log(exp(-1.5 * x) + 1) #default drift auto effect transform

#raw pop sd transform
sdscale <- 1 #default
rawsdtform <- function(x) exp(x * 2 -2) * sdscale #default

#sd approximation function
sdapprox <- function(means,sds,tform) {
  for(i in 1:length(means)){
    sds[i] <- ((tform(means[i]+sds[i]*3) - tform(means[i]-sds[i]*3))/6 + 
      (tform(means[i]+sds[i]) - tform(means[i]-sds[i]))/2) /2
  }
  return(sds)
}

#raw population mean parameters
rawpopmeans_prior <- rnorm(n, 0, 1) #prior distribution for rawpopmeans
rawpopmeans_sample <- -.3 #hypothetical sample
sdscale <- 1 #default

#population mean parameters after parameter specific transform
popmeans_prior <- tform(rawpopmeans_prior)
popmeans_sample <- tform(rawpopmeans_sample)

#plot pop means
plot(density(rawpopmeans_prior), ylim=c(0,1), xlim=c(-5,2),
  xlab='Parameter value', main='Population means')
points(density(popmeans_prior, bw=bw),col=2,type='l')
segments(y0=0,y1=.5,x0=c(rawpopmeans_sample,popmeans_sample),lty=3,col=1:2) 
legend('topleft',c('Raw pop. mean prior', 'Pop. mean prior', 
  'Raw pop. mean sample', 'Pop. mean sample'),lty=c(1,1,3,3), col=1:2, bty='n')

#population standard deviation parameters
rawpopsd_prior <- rawsdtform(rnorm(n, 0, 1)) #raw population sd prior

popsd_prior <- sdapprox(rawpopmeans_prior,rawpopsd_prior,tform)

#sample population standard deviation posterior
rawpopsd_sample <- rawsdtform(.9) #hypothetical sample
popsd_sample <- sdapprox(means=rawpopmeans_sample,  #transform sample to actual pop sd
  sds=rawpopsd_sample,tform=tform)

#plot pop sd
plot(density(rawpopsd_prior,from=-.2,to=10,na.rm=TRUE, bw=bw), xlab='Parameter value', 
  xlim=c(-.1,3), ylim=c(0,2), main='Population sd')
points(density(popsd_prior,from=-.2,to=10,na.rm=TRUE, bw=bw),type='l', col=2)
segments(y0=0,y1=1,x0=c(rawpopsd_sample, popsd_sample), col=1:2,lty=3)
legend('topright',c('Raw pop. sd prior','Pop. sd prior', 
  'Raw pop. sd sample','Pop. sd sample'), col=1:2, lty=c(1,1,3,3),bty='n')

#individual level parameters

#marginal individual level parameters (given all possible values for mean and sd)
rawindparams_margprior <- rawpopmeans_prior + rawpopsd_prior * rnorm(n, 0, 1) 
indparams_margprior <- tform(rawindparams_margprior)

plot(density(rawindparams_margprior,from=-10,to=10,bw=bw), xlab='Parameter value', 
  xlim=c(-5,2), ylim=c(0,1), main='Marginal dist. individual parameters')
points(density(indparams_margprior,from=-10,to=.2,bw=bw),type='l',col=2)
legend('topleft',c('Raw individual parameters prior','Individual parameters prior'), 
  col=1:2,lty=1,bty='n')


#conditional individual level parameters (given sampled values for mean and sd)
rawindparams_condprior<- rawpopmeans_sample + rawpopsd_sample * rnorm(n,0,1) 
rawindparams_condsample<- rawpopmeans_sample + rawpopsd_sample * rnorm(nsubjects,0,1) 
indparams_condprior<- tform(rawindparams_condprior)
indparams_condsample<- tform(rawindparams_condsample)


plot(density(rawindparams_condprior), xlab='Parameter value', xlim=c(-5,2),
  ylim=c(0,1), main='Conditional dist. individual parameters')
points(density(indparams_condprior),type='l',col=2)
segments(y0=0,y1=.5,x0=c(rawindparams_condsample, indparams_condsample),
  col=rep(1:2,each=nsubjects),lty=3, lwd=2)
legend('topleft',c('Raw ind. pars. prior','Ind. pars. prior', 
  'Raw ind. pars. samples','Ind. pars. samples'), col=1:2, lty=c(1,1,3,3),bty='n')


