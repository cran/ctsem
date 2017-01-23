if(1==0){
out=array(NA,dim = c(19,9,5))
for(iteri in 1:dim(out)[3]){
  for(nsubjects in c(200)){
    for(tpoints in c(10)){
      
      nlatent=2
      nmanifest=2
      
      genm=ctModel(Tpoints=tpoints,
        n.latent=nlatent, n.manifest=nmanifest, 
        manifestNames=c('man1','man2'),
        LAMBDA=matrix(c(1, 0,0,1), nrow=nmanifest, ncol=nlatent),
        DRIFT=matrix(c(-.4, .1, 0, -0.2), byrow=TRUE, nrow=nlatent, ncol=nlatent),
        DIFFUSION=t(chol(matrix(c(2, 0, 0, 2), byrow=TRUE, nrow=nlatent, ncol=nlatent))),
        MANIFESTVAR=t(chol(matrix(c(.4, 0,0,.4), nrow=nmanifest, ncol=nmanifest))),
        CINT=matrix(c(0, 0), nrow=nlatent, ncol=1),
        MANIFESTTRAITVAR=t(chol(matrix(c(6,3,3,6),nmanifest,nmanifest))),
        MANIFESTMEANS=matrix(c(0,0), nrow=nmanifest, ncol=1))
      
      dat=ctGenerate(ctmodelobj=genm, n.subjects=nsubjects, burnin=30, dtmean=.9, 
        logdtsd=0,simultdpredeffect=TRUE,wide=TRUE)
      # plot(dat[1,1:tpoints])
      ctIndplot(datawide = dat,n.subjects = 3,n.manifest = nmanifest,Tpoints = tpoints)
      
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
      ftime=tb-ta
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
      # sm$pars$indvarying[sm$pars$matrix =='DRIFT']=TRUE
      # sm$pars$transform[sm$pars$matrix =='DRIFT' & sm$pars$row == sm$pars$col]='param*5-1'
      
      ta=proc.time()[3]
      sfit=ctStanFit(datalong=ldat, ctstanmodel=sm,iter=300,chains=2,esthyper = TRUE,optimize=T)
      tb=proc.time()[3]
      btime=tb-ta
      sum2=summary(sfit)$popmeans
      rownames(sum2)=gsub('hmean_','',rownames(sum2))
      sum2=sum2[match(rownames(est),rownames(sum2)),]
      sum2[c(17,19),]=summary(sfit)$popsd
      
      tmpdir=tempdir()
      tmpdir=gsub('\\','/',tmpdir,fixed=T)
      system(paste0("rm ",tmpdir,'/*.*'))
      
      sm=ctStanModel(ctmodelobj = fitm,type = 'stanct')
      sm$pars$indvarying=TRUE
      # sm$pars$indvarying[sm$pars$matrix =='DRIFT']=TRUE
      # sm$pars$transform[sm$pars$matrix =='DRIFT' & sm$pars$row == sm$pars$col]='param*5-1'
      
      ta=proc.time()[3]
      sfit=ctStanFit(datalong=ldat, ctstanmodel=sm,iter=300,chains=2,esthyper = TRUE,optimize=F)
      tb=proc.time()[3]
      bstime=tb-ta
      sum3=summary(sfit)$popmeans
      rownames(sum3)=gsub('hmean_','',rownames(sum3))
      sum3=sum3[match(rownames(est),rownames(sum3)),]
      sum3[c(17,19),]=summary(sfit)$popsd
      
      tmpdir=tempdir()
      tmpdir=gsub('\\','/',tmpdir,fixed=T)
      system(paste0("rm ",tmpdir,'/*.*'))
      
      ta=proc.time()[3]
      sm$pars$indvarying[sm$pars$matrix =='DRIFT']=TRUE
      sfit=ctStanFit(datalong=ldat, ctstanmodel=sm,iter=300,chains=2,esthyper = T,optimize=F)
      tb=proc.time()[3]
      bsnhtime=tb-ta
      sum4=summary(sfit)$popmeans
      rownames(sum4)=gsub('hmean_','',rownames(sum4))
      sum4=sum4[match(rownames(est),rownames(sum4)),]
      sum4[c(17,19),]=summary(sfit)$popsd[c(5,6),]
      
      tmpdir=tempdir()
      tmpdir=gsub('\\','/',tmpdir,fixed=T)
      system(paste0("rm ",tmpdir,'/*.*'))
      
      outmat=cbind(est[,1],sum2[,1],outmatbak[,3],sum4[,3],
        est[,2],sum2[,2],outmatbak[,6],sum4[,4],
        est[,3],sum2[,3],outmatbak[,9],sum4[,5])
      
      colnames(outmat)=1:ncol(outmat)
      colnames(outmat)[1:12]=c(paste0(rep(c('lb ','est ','ub '),each=4), 
        c('freq', 'optim', 'hmc', 'hmc_nohyper')))
      
      outmat=round(outmat,3)
      
      print(outmat)
      
      print(paste0('frequentist seconds = ',ftime,'   Bayes optim seconds = ', btime, 
        '   Bayes sample seconds = ',bstime, '    Sample no hyper seconds = ', bsnhtime))
      out[,,iteri]=outmat
      
    }
  }
}
}
