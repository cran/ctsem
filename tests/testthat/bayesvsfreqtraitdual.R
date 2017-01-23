if(1==0){
out=array(NA,dim = c(19,15,5))
for(iteri in 1:dim(out)[3]){
  for(nsubjects in c(5,20)){
    for(tpoints in c(6,25)){
      
      nlatent=2
      nmanifest=2
      
      genm=ctModel(Tpoints=tpoints,
        n.latent=nlatent, n.manifest=nmanifest, 
        manifestNames=c('man1','man2'),
        LAMBDA=matrix(c(1, 0,0,1), nrow=nmanifest, ncol=nlatent),
        DRIFT=matrix(c(-.5, .1, 0, -0.2), byrow=TRUE, nrow=nlatent, ncol=nlatent),
        DIFFUSION=matrix(c(3, -2, -2, 4), byrow=TRUE, nrow=nlatent, ncol=nlatent),
        MANIFESTVAR=matrix(c(1.5, 0,0,.7), nrow=nmanifest, ncol=nmanifest),
        CINT=matrix(c(0, 0), nrow=nlatent, ncol=1),
        MANIFESTTRAITVAR=matrix(c(2,.5,0,2),nmanifest,nmanifest),
        MANIFESTMEANS=matrix(c(-5,2), nrow=nmanifest, ncol=1))
      
      dat=ctGenerate(ctmodelobj=genm, n.subjects=nsubjects, burnin=0, dtmean=.9, 
        logdtsd=0,simultdpredeffect=TRUE,wide=TRUE)
      # plot(dat[1,1:tpoints])
      ctIndplot(datawide = dat,n.subjects = 1,n.manifest = nmanifest,Tpoints = tpoints)
      
      fitm= ctModel(Tpoints=tpoints, manifestNames=c('man1','man2'),
        MANIFESTTRAITVAR='auto',
        n.latent=nlatent, n.manifest=nmanifest,LAMBDA=diag(nlatent)
      )
      
      dat[1,'man1_T0']=NA
      dat[2,paste0(fitm$manifestNames,'_T',rep(0:(tpoints-1),each=2))]=NA
      
      ta=proc.time()[3]
      fit=ctFit(dat, fitm,verbose=0,nofit=F,retryattempts=1) #,objective='Kalman')
      summary(fit)$ctparameters
      tb=proc.time()[3]
      # fit=ctCI(ctfitobj = fit,confidenceintervals = pars[1,pari])
      # tc=proc.time()[3]
      
      se1=summary(fit,verbose=TRUE)$omxsummary$parameters
      varpars=se1[,1][se1[,'row']==se1[,'col']][grep(
        'diffusion|T0var|manifestvar|manifesttrait',se1[,1][se1[,'row']==se1[,'col']])]
      est=cbind(se1[,5]-1.96*se1[,6],se1[,5],se1[,5]+1.96*se1[,6])
      est[se1[,1] %in% varpars,] = exp( est[se1[,1] %in% varpars,])
      rownames(est)=se1[,1]
      colnames(est)=c('lower','estimate','upper')
      est=round(est,3)
      est[est>999]<-999
      # est1=fit$mxobj$output$estimate[names(fit$mxobj$output$estimate)==pars[1,pari]]
      # se=c(est1 - se1*2, est1,est1+2*se1,NA,fit$mxobj$output$status$code,tb-ta)
      
      # ci=c(fit$mxobj$output$confidenceIntervals,NA,fit$mxobj$output$status$code,tc-ta)
      # if(varpar) ci = exp(ci)
      
      
      # expect_equal(as.numeric(summary(fit)$ctparameters[1]),.4,tolerance=.02)
      
      ldat=ctWideToLong(datawide = dat,Tpoints = tpoints,n.manifest = nmanifest,n.TDpred = 0,n.TIpred = 0,
        manifestNames=fitm$manifestNames)
      ldat=ctDeintervalise(ldat)
      
      sm=ctStanModel(ctmodelobj = fitm,type = 'stanct')
      sm$pars$indvarying[sm$pars$matrix !='MANIFESTMEANS']=FALSE
      
      ta=proc.time()[3]
      sfit=ctStanFit(datalong=ldat, ctstanmodel=sm,iter=300,chains=2,esthyper = TRUE,optimize=F)
      tb=proc.time()[3]
      sum2=summary(sfit)$popmeans[,c('2.5%','50%','97.5%','n_eff','Rhat')]
      rownames(sum2)=gsub('hmean_','',rownames(sum2))
      sum2=sum2[match(rownames(est),rownames(sum2)),]
      sum2[c(17,19),]=summary(sfit)$popsd[,c('2.5%','50%','97.5%','n_eff','Rhat')]
      
      # khdi=c(summary(sfit)$popmeans[rownames(summary(sfit)$popmeans)==
      #     paste0('hmean_',pars[1,pari]),c('2.5%','50%','97.5%','n_eff','Rhat')],tb-ta)
      
      tmpdir=tempdir()
      tmpdir=gsub('\\','/',tmpdir,fixed=T)
      system(paste0("rm ",tmpdir,'/*.*'))
      
      ta=proc.time()[3]
      sfit3=ctStanFit(datalong=ldat, ctstanmodel=sm,iter=300,chains=2,esthyper = F,optimize=T,kalman=TRUE)
      sum3=summary(sfit3)$popmeans[,c('2.5%','50%','97.5%','n_eff','Rhat')]
      rownames(sum3)=gsub('hmean_','',rownames(sum3))
      sum3=sum3[match(rownames(est),rownames(sum3)),]
      sum3[c(17,19),]=summary(sfit3)$popsd[,c('2.5%','50%','97.5%','n_eff','Rhat')]
      tb=proc.time()[3]
      # shdi=c(summary(sfit)$popmeans[rownames(summary(sfit)$popmeans)==
      #     paste0('hmean_',pars[1,pari]),c('2.5%','50%','97.5%','n_eff','Rhat')],tb-ta)
      
      tmpdir=tempdir()
      tmpdir=gsub('\\','/',tmpdir,fixed=T)
      system(paste0("rm ",tmpdir,'/*.*'))
      
      # row = cbind(as.numeric(pars[2,pari]),round(rbind(se, ci,khdi,shdi),3))
      # rownames(row)[1]=paste(pars[1,pari],', N=',nsubjects,'T=',tpoints,', ','stdError')
      # rownames(row)[2]=paste(pars[1,pari],', N=',nsubjects,'T=',tpoints,', ','ci')
      # rownames(row)[3]=paste(pars[1,pari],', N=',nsubjects,'T=',tpoints,', ','esthyper hdi')
      # rownames(row)[4]=paste(pars[1,pari],', N=',nsubjects,'T=',tpoints,', ','!esthyper hdi')
      
      # expect_equal(as.numeric(ci[2]),as.numeric(hdi[2]),tolerance=.03)
      
      # out<-rbind(out,row)
      # colnames(out)=c('true','lower','estimate','upper','n_eff','Code/Rhat','seconds')
      # print(out)
      
      outmat=cbind(est[,1],sum2[,1],sum3[,1],
        est[,2],sum2[,2],sum3[,2],
        est[,3],sum2[,3],sum3[,3],
        NA,sum2[,4],sum3[,4],
        NA,sum2[,5],sum3[,5])
      colnames(outmat)=1:ncol(outmat)
      colnames(outmat)[1:9]=c(paste0(rep(c('lb ','est ','ub '),each=3), 
        c('StdError', 'esthyper hdi', '!esthyper hdi')))
      colnames(outmat)[10:15]=c(paste0(rep(c('n_eff ','Rhat '),each=3), 
        c('StdError', 'esthyper hdi', '!esthyper hdi')))
      
      print(outmat)
      out[,,iteri]=outmat
      
    }
  }
}
}
