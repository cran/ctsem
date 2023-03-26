if(1==99){
  
  checkm<-ctModel(
    type='stanct',
    n.latent=2,n.TDpred=1,n.TIpred=1,n.manifest=2,
    MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
    MANIFESTMEANS=0,
    DRIFT=c('dr1','dr12','dr21||||TI1','dr22'),
    DIFFUSION=c('diff11',0,'diff21','diff22||||TI1'),
    CINT=matrix(c('cint1||||TI1','cint2||||TI1'),ncol=1),
    LAMBDA=diag(2),tipredDefault=FALSE)  
  
  ctstantestfit<-ctStanFit(ctsem::ctstantestdat,checkm,cores=1,control=list(max_treedepth=4),
    optimize = F,optimcontrol=list(finishsamples=20,stochastic=F),nopriors=FALSE)
  
  s=ctstantestfit$standata
  pars <- ctstantestfit$stanfit$rawest
  
  Jkalman <- function(s, pars){
    
    ms <- data.frame(s$matsetup)
    mv <- data.frame(s$matvalues)
    
    #preallocate
    II <- diag(1,s$nlatentpop)
    
    yprior <- rep(NA, s$nmanifest)
    ycov <- matrix(NA, s$nmanifest,s$nmanifest)
    ychol <- matrix(NA, s$nmanifest,s$nmanifest)
    state <- rep(NA,s$nlatentpop)
    statecov <- matrix(NA, s$nlatentpop,s$nlatentpop)
    
    
    
    T0VAR <- T0cov <- JA <- matrix(NA, s$nlatentpop,s$nlatentpop)
    A <- dtA <- matrix(NA, s$nlatent,s$nlatent)
    T0MEANS <- CINT <- dtb <- matrix(NA, s$nlatent)
    
    TDPREDEFFECT <- matrix(NA, s$nlatent, s$ntdpred)
    JTD <- matrix(NA, s$nlatentpop, s$ntdpred)
    
    G <- matrix(NA, s$nlatent,s$nlatent)
    Q <- matrix(NA, s$nlatent,s$nlatent)
    dtQ <- matrix(NA, s$nlatent,s$nlatent)
    Qinf <- matrix(NA, s$nlatent,s$nlatent)
    
    di <- 1:s$nlatent
    
    LAMBDA <- matrix(NA, s$nmanifest,s$nlatent)
    JLAMBDA <- matrix(NA, s$nmanifest, s$nlatentpop)
    MANIFESTMEANS <- matrix(NA, s$nmanifest)
    MANIFESTVAR <- Mcov <- matrix(NA,s$nmanifest,s$nmanifest)
    
    
    K <- matrix(NA, s$nlatent,s$nmanifest)
    ll <- rep(NA, s$ndatapoints)
    
    t0step <- function(){
      state[di] <- T0MEANS
      statecov <- T0cov
    }
    
    dynamicstep <- function(){
      state <<- dtA %*% state[di] + dtb
      statecov <<- JdtA %*% statecov %*% t(JdtA) 
      statecov[di,di] <<- statecov[di,di] + dtQ
    }
    
    tdpredstep <- function(){
      state[di] <<- state[di] + TDPREDEFFECT %*% tdpreds[i,]
      statecov <<- JTD %*% statecov %*% t(JTD)
    }
    
    stateupdate <- function(){
      yprior[oi,] <<- LAMBDA[oi, ] %*% state[di] + MANIFESTMEANS[oi]
      ycov[oi,oi] <<- JLAMBDA[oi, ] %*% statecov %*% t(JLAMBDA[oi, oi]) + Mcov[oi,oi]
      K[,oi] <<- ycov[oi,oi] %*% t(JLAMBDA[oi,]) %*% solve(ycov[oi,oi])
      state <<- state + K[,oi] %*% y[i,]
      statecov <<- (II - K %*% JLAMBDA) %*% statecov
    }
    
    loglik <- function(){
      return( dnorm(di = y[i,oi] - yprior[oi], log=TRUE) + log(trace(ychol[oi,oi])))
    }
    
    
    
    dtAfunc <- function(dt) dtA<<-expm::expm(A * dt)
    dtJAfunc <- function(dt) dtJA<<-expm::expm(JA * dt)
    dtbfunc <- function() dtb <<- Ainv %*% (dtA-II[di,di]) %*% CINT
    
    Qinffunc <- function(){
      Ahatch<<-A %di% II[di,di] +  II[di,di] %di% A
      Qinf<<-matrix(-solve(Ahatch , Q), nrow=nrow(A))
    }
    
    dtQfunc <- function(Qinf, dtJA) dtQ <<- Qinf - (dtJA %*% Qinf %*% t(dtJA ))
    
    fillSysMats<-function(when){
      for(ri in 1:nrow(ms)){
        if(ms$when[ri]==when){
          if(when > 0 || ms$indvarying[ri]>0 || ms$tipred[ri] >0 || si==1){
            
            
          }
        }
      }
    }
    
    tformState<-function(){
      
    }
    
    
    si <- t0check <- 0
    for(i in 1:ndatapoints){
      
      if(s$subject[i] != si) t0check <- 0 else t0check <- t0check + 1
      si <- s$subject[i]
      
      if(t0check > 0) dt <- s$time[i]-prevtime
      prevtime <- s$time[i]
      
      #t0 setup
      if(t0check==0){
        fillSysMats(when=0)
        fillSysMats(when=1)
        t0step()
      }
      
      if(t0check > 0){
        fillSysMats(when=2)
        dtAfunc()
        dtbfunc()
        Qinffunc()
        dtQfunc()
        dynamicstep()
      }
      
      if(s$ntdpred > 0){
        fillSysMats(when=3)
        tdpredstep()
      }
      
      oi <- s$which
      if(length(oi) > 0){
        stateupdate()
        ll[i]<-loglik()
      }
      
    }
    
    return(sum(ll,na.rm=TRUE))
  }
  
  
  
  
  
}
