if(1==0){
  library(ctsem)
  
  sunspot<-sunspot.year
  sunspot<-sunspot[50: (length(sunspot)-(1988-1924))]
  id <- 1
  time <- 1749:1924
  datalong <- cbind(id, time, sunspot)
  
  TT=176
  
  model1 <- ctModel(type='omx', n.latent = 2, n.manifest = 1, Tpoints = TT, 
    LAMBDA = matrix(c(1, .6), nrow = 1, ncol = 2,byrow=TRUE),
    manifestNames = 'sunspot',
    DRIFT = matrix(c(
      0,  1,
      -.3,   -.35
    ),nrow =2, ncol = 2,byrow=TRUE),
    MANIFESTMEANS = matrix(c(46),nrow=1,ncol=1),
    MANIFESTVAR = matrix(c(
      0), nrow=1, ncol=1),
    CINT =matrix(c(0,0),nrow=2,ncol=1),
    T0MEANS=matrix(c(0,0),nrow=2,ncol=1),
    T0VAR=diag(10,2),
    PARS=matrix('p'),
    DIFFUSION=matrix(c(0,0,
      0,16),ncol=2,nrow=2,byrow=TRUE))
  
  
  
  summary(fit1$mxobj)
  sqrt(fit1$mxobj$DIFFUSION$values[2,2]) 
  
  
  sm <- ctStanModel(model1)
  
  ssfit<-ctStanFit(datalong=datalong, ctstanmodel=sm, forcerecompile = T,
    verbose=0,optimize=T,#optimcontrol=list(estonly=T),
    # nlcontrol=list(nldynamics=F, nlmeasurement=F,ukffull=1,ukfspread=1e-1),
    chains=1,priors=TRUE)
  
}

