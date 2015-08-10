## ----setup, include = FALSE, cache = FALSE, echo = FALSE----------------------
library('ctsem')
library(knitr)
render_sweave()
set.seed(22)
opts_chunk$set(fig.path = 'figures/plots-', warning = FALSE, fig.align = 'center', width.cutoff = 80, fig.show = 'hold', eval = TRUE, echo = TRUE, message = FALSE, background = "white", prompt = TRUE, highlight = FALSE, comment = NA, tidy = FALSE, out.truncate = 80)
options(replace.assign = TRUE, width = 80, prompt = "R> ", scipen = 12, digits = 3)



# setwd('C:\\Users\\driver\\Dropbox\\MPIB\\CT-SEM\\manual') #set this working directory!
Sys.setenv(TEXINPUTS = getwd(),
  BIBINPUTS = getwd(),
  BSTINPUTS = getwd())

## ----install, echo = T, eval = F----------------------------------------------
#  install.packages("ctsem", repos = "http://r-forge.r-project.org", type = 'source')
#  library('ctsem')

## ----wideformat, echo = FALSE-------------------------------------------------
data('datastructure')
datastructure

## ----longformat, include = TRUE, cache = FALSE, echo = FALSE, results = 'markup'----
data('longexample')
head(longexample, 7)

## ----longformatconversion, include = TRUE, cache = FALSE, echo = TRUE, results = 'markup'----
data('longexample')
wideexample <- ctLongToWide(datalong = longexample, id = "subject", 
  time = "Time", manifestNames = c("Y1", "Y2", "Y3"), 
  TDpredNames = "TD1", TIpredNames = c("TI1", "TI2"))
wide <- ctIntervalise(datawide = wideexample, Tpoints = 3, n.manifest = 3, 
  n.TDpred = 1, n.TIpred = 2, manifestNames = c("Y1", "Y2", "Y3"), 
  TDpredNames = "TD1", TIpredNames = c("TI1", "TI2") )

## ----simplemodel, include = TRUE, echo = TRUE, results = 'hide'---------------
examplemodel <- ctModel(n.latent = 2, n.manifest = 2, Tpoints = 3, 
  LAMBDA = diag(2))

## ----example1ctfit, include = TRUE, cache = TRUE, echo = TRUE, results = 'hide'----
data('ctExample1')
example1model <- ctModel(n.latent = 2, n.manifest = 2, Tpoints = 6, 
  manifestNames = c('LeisureTime', 'Happiness'), 
  latentNames = c('LeisureTime', 'Happiness'), LAMBDA = diag(2))
example1fit <- ctFit(datawide = ctExample1, ctmodelobj = example1model)
summary(example1fit)

## ----example1ctfittable, include = TRUE, echo = TRUE--------------------------
summary(example1fit)['discreteDRIFTstd']

## ----example1testing, cache = TRUE, echo = TRUE-------------------------------
testmodel <- example1model
testmodel$DRIFT[1, 2] <- 0
testfit <- ctFit(datawide = ctExample1, ctmodelobj = testmodel)

## ----mxcompare----------------------------------------------------------------
mxCompare(example1fit$mxobj, testfit$mxobj)

## ----confidenceintervals, cache = TRUE, echo = 1------------------------------
example1cifit <- ctFit(datawide = ctExample1, ctmodelobj = example1model, 
  confidenceintervals = 'DRIFT')
example1cifit$mxobj$output$confidenceIntervals

## ----example2fit, cache = TRUE------------------------------------------------
data('ctExample1')
traitmodel <- ctModel(n.manifest = 2, n.latent = 2, Tpoints = 6, 
  LAMBDA = diag(2), manifestNames = c('LeisureTime', 'Happiness'), 
  latentNames = c('LeisureTime', 'Happiness'), TRAITVAR = "auto")
traitfit <- ctFit(datawide = ctExample1, ctmodelobj = traitmodel)

## ----traitparamplot, include = TRUE, cache = FALSE, echo = FALSE, results = 'hide', fig.height = 4----
par(mfrow = c(2, 2))
par(mar = c(2, 4, 2, 2))
plot(example1fit, wait = FALSE, max.time = 20, mean = FALSE, withinVariance = FALSE, betweenVariance = FALSE)
plot(traitfit, wait = FALSE, max.time = 20, mean = FALSE, withinVariance = FALSE, betweenVariance = FALSE)

## ----example1TIpred, include = TRUE, cache = TRUE, echo = TRUE, results = 'hide'----
data('ctExample1TIpred')
tipredmodel <- ctModel(n.manifest = 2, n.latent = 2, n.TIpred = 1,
  manifestNames = c('LeisureTime', 'Happiness'),
  latentNames = c('LeisureTime', 'Happiness'),
  TIpredNames = 'NumFriends',
 Tpoints = 6, LAMBDA = diag(2), TRAITVAR = "auto")
tipredfit <- ctFit(datawide = ctExample1TIpred, ctmodelobj = tipredmodel)

summary(tipredfit)['TIPREDEFFECT']
summary(tipredfit)['discreteTIPREDEFFECT']
summary(tipredfit)['asymTIPREDEFFECT']
summary(tipredfit)['addedTIPREDVAR']

## ----example1TIpredestimates, echo = FALSE, out.width = '4cm', out.height = '4cm'----
summary(tipredfit)['TIPREDEFFECT']
cat('\n')
summary(tipredfit)['discreteTIPREDEFFECT']

## ----example1TIpredestimates2, echo = FALSE, out.width = '4cm', out.height = '4cm'----
summary(tipredfit)['asymTIPREDEFFECT']
cat('\n')
summary(tipredfit)['addedTIPREDVAR']

## ----tdpreddemo, echo = FALSE, eval = TRUE, fig.height = 3--------------------
Tpoints = 20
testm <- ctModel(Tpoints = Tpoints, n.latent = 1, n.TDpred = 1, n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.1, 1), TDPREDEFFECT = diag(1, 1), TDPREDVAR = diag(0, Tpoints-1), CINT = diag(4, 1), T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(100, 1), TDPREDMEANS = matrix(c(0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ncol = 1, nrow = (Tpoints-1)))
testd <- ctGenerate(testm, n.subjects = 100, burnin = 300)

testm <- ctModel(Tpoints = Tpoints, n.latent = 1, n.TDpred = 1, n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.0, 1), TDPREDEFFECT = diag(1, 1), TDPREDVAR = diag(0, Tpoints-1), CINT = diag(4, 1), T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(1, 1), TDPREDMEANS = matrix(c(0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ncol = 1, nrow = (Tpoints-1)))
testdpure <- ctGenerate(testm, n.subjects = 1, burnin = 300)

par(mfrow = c(1, 2), cex = .7)
plot(0:(Tpoints-1), testdpure[, 1:Tpoints], type = 'l', ylim = c(min(testd[, 1:Tpoints]), max(testd[, 1:Tpoints])),
  ylab = 'Dependent variable', xlab = 'Time', lwd = 3, main = 'Impulse predictor')
for(i in 1:5){
  points(0:(Tpoints-1), testd[i, 1:Tpoints], col = 1+i, type = 'b')
}

Tpoints = 20
testm <- ctModel(Tpoints = Tpoints, n.latent = 1, n.TDpred = 1, n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.15, 1), TDPREDEFFECT = diag(1.6, 1), TDPREDVAR = diag(0, Tpoints-1), CINT = diag(4, 1), T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(100, 1), TDPREDMEANS = matrix(c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), ncol = 1, nrow = (19)))
testd <- ctGenerate(testm, n.subjects = 100, burnin = 300)

testm <- ctModel(Tpoints = Tpoints, n.latent = 1, n.TDpred = 1, n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.0, 1), TDPREDEFFECT = diag(1.6, 1), TDPREDVAR = diag(0, Tpoints-1), CINT = diag(4, 1), T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(1, 1), TDPREDMEANS = matrix(c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), ncol = 1, nrow = (19)))
testdpure <- ctGenerate(testm, n.subjects = 1, burnin = 300)

plot(0:(Tpoints-1), testdpure[, 1:Tpoints], type = 'l', ylim = c(min(testd[, 1:Tpoints]), max(testd[, 1:Tpoints])), ylab = 'Dependent variable', xlab = 'Time', lwd = 3, main = 'Level predictor')
for(i in 1:5){
  points(0:(Tpoints-1), testd[i, 1:Tpoints], col = 1+i, type = 'b')
}

## ----example2TDpred, include = TRUE, cache = TRUE, echo = TRUE, results = 'hide', fig.height = 4, fig.width = 4, fig.align = 'center'----
data('ctExample2')
tdpredmodel <- ctModel(n.manifest = 2, n.latent = 2, n.TDpred = 1, 
  Tpoints = 8, manifestNames = c('LeisureTime', 'Happiness'), 
  TDpredNames = 'MoneyInt', latentNames = c('LeisureTime', 'Happiness'),
  T0TDPREDCOV = matrix(0, nrow = 2, ncol=7),
  TRAITTDPREDCOV = matrix(0, nrow = 2, ncol=7), 
  LAMBDA = diag(2), TRAITVAR = "auto")
tdpredfit <- ctFit(datawide = ctExample2, ctmodelobj = tdpredmodel)

summary(tdpredfit)['TDPREDEFFECT']
summary(tdpredfit)['discreteTDPREDEFFECT']

## ----example2TDpredestimates, echo = FALSE, out.width = '3cm'-----------------
summary(tdpredfit)['TDPREDEFFECT']

## ----example2TDpredestimates3, echo = FALSE, out.width = '3cm'----------------
summary(tdpredfit)['discreteTDPREDEFFECT']

## ----example2TDpredlevel, include = TRUE, eval = TRUE, echo = TRUE, cache = TRUE----
data('ctExample2')
tdpredmodel <- ctModel(n.manifest = 2, n.latent = 3, n.TDpred = 1, 
  Tpoints = 8, manifestNames = c('LeisureTime', 'Happiness'), 
  TDpredNames = 'MoneyInt', 
  latentNames = c('LeisureTime', 'Happiness', 'MoneyIntLatent'),
  T0TDPREDCOV = matrix(0, nrow = 3, ncol = 7),
  TRAITTDPREDCOV = matrix(0, nrow = 3, ncol = 7), 
  LAMBDA = matrix(c(1,0, 0,1, 0,0), ncol = 3), TRAITVAR = "auto")

tdpredmodel$TRAITVAR[3, ] <- 0
tdpredmodel$TRAITVAR[, 3] <- 0
tdpredmodel$DIFFUSION[, 3] <- 0
tdpredmodel$DIFFUSION[3, ] <- 0
tdpredmodel$T0VAR[3, ] <- 0
tdpredmodel$T0VAR[, 3] <- 0
tdpredmodel$CINT[3] <- 0
tdpredmodel$T0MEANS[3] <- 0
tdpredmodel$TDPREDEFFECT[3, ] <- 1
tdpredmodel$DRIFT[3, ] <- 0

tdpredfit <- ctFit(datawide = ctExample2, ctmodelobj = tdpredmodel)

summary(tdpredfit)['DRIFT']
summary(tdpredfit, timeInterval = 20)['discreteTDPREDEFFECT']

## ----timeseries, cache = TRUE, echo = TRUE------------------------------------
data('ctExample3')
model <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 100, 
  LAMBDA = matrix(c(1, 'lambda2', 'lambda3'), nrow = 3, ncol = 1), 
  MANIFESTMEANS = matrix(c(0, 'manifestmean2', 'manifestmean3'), nrow = 3, 
    ncol = 1))
fit <- ctFit(data = ctExample3, ctmodelobj = model, objective = 'Kalman', 
  stationary = c('T0VAR'))

## ----multigroup, cache = TRUE, echo = TRUE------------------------------------
data('ctExample4')

basemodel <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 20,
  LAMBDA = matrix(c(1, 'lambda2', 'lambda3'), nrow = 3, ncol = 1),
  TRAITVAR='auto', MANIFESTMEANS = matrix(c(0, 'manifestmean2', 
    'manifestmean3'), nrow = 3, ncol = 1))

freemodel <- basemodel
freemodel$LAMBDA[3, 1] <- 'groupfree'
groups <- paste0('g', rep(1:2, each = 10), '_')

multif <- ctMultigroupFit(datawide = ctExample4, groupings = groups,
  ctmodelobj = basemodel, freemodel = freemodel)

## ----multigroupOutput, echo = FALSE-------------------------------------------
multif$output$estimate[grep('lambda3', names(multif$output$estimate))]

## ----dynamicresidualmodel, echo = TRUE, results='hide',fig.keep='none'--------
testm <- ctModel(Tpoints = 200, n.latent = 2, n.manifest = 1, 
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  DIFFUSION = matrix(c(0, 0, 0, 1), 2),
  MANIFESTVAR = diag(.4,1),
  DRIFT = matrix(c(0, -.1, 1, -.3), nrow = 2),   
  CINT = matrix(c(1, 0), nrow = 2))

data<-ctGenerate(testm,n.subjects=1,burnin=200,dT=1)

ctIndplot(data,n.subjects=1,n.manifest=1,Tpoints=200)

model <- ctModel(Tpoints = 200, n.latent = 2, n.manifest = 1, 
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  DIFFUSION = matrix(c(0, 0, 0, 'diffusion'), 2),
  DRIFT = matrix(c(0, 'regulation', 1, "diffusionAR"), nrow = 2),   
  CINT = matrix(c("processCINT", 0), nrow = 2))

fit<-ctFit(data,model,stationary=c('T0MEANS','T0VAR'))

## ----osscilating, cache = TRUE, echo = TRUE, eval = TRUE, include = TRUE------
data('Oscillating')

inits <- c(-38, -.5, 1, 1, .1, 1, 0, .9)
names(inits) <- c('cross','auto', 'diffusion22',
  'T0var11', 'T0var21', 'T0var22','m1', 'm2')

oscillatingm <- ctModel(n.latent = 2, n.manifest = 1, Tpoints = 11, 
  MANIFESTVAR = matrix(c(0), nrow = 1, ncol = 1), 
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  T0VAR = matrix(c('T0var11', 'T0var21', 0, 'T0var22'), nrow = 2, ncol = 2),
  DRIFT = matrix(c(0, "crosseffect", 1, "autoeffect"), nrow = 2, ncol = 2), 
  CINT = matrix(0, ncol = 1, nrow = 2),
  DIFFUSION = matrix(c(0, 0, 0, "diffusion"), nrow = 2, ncol = 2),
  startValues = inits)

oscillatingf <- ctFit(Oscillating, oscillatingm)

## ----inits, echo = TRUE-------------------------------------------------------
omxInits <- omxGetParameters(example1fit$mxobj)

fitWithInits <- ctFit(data = ctExample1, ctmodelobj = example1model, 
  omxStartValues = omxInits)

