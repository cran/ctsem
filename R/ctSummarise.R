if(1==99){
ctSummarise<-function(sf,name='ctSummary',cores=2, times=seq(0,10,.1),
  latents=1:sf$ctstanmodelbase$n.latent, manifests=1:sf$ctstanmodelbase$n.manifest){

if(class(sf)=='ctStanFit'){ #avoid overwriting plots if error!
  library(data.table)

  ydat=ctsem:::standatatolong(sf$standata,origstructure = TRUE,ctm=sf$ctstanmodelbase)

  sm <- sf$ctstanmodelbase 
  s=summary(sf)
  if(is.null(sf$generated)) sf <- ctStanGenerateFromFit(fit = sf,nsamples = 100,fullposterior = FALSE,cores = cores)
  cp <- ctStanContinuousPars(sf)
  
  n<-sapply(unique(ydat[[sm$subjectIDname]]),function(x) sum(ydat[[sm$subjectIDname]] %in% x))
  whichsubfull <- unique(ydat[[sm$subjectIDname]])[order(n,decreasing = TRUE)][1:(min(30,length(n)))]
  # whichsubfull <- match(whichsubfull,sf$setup$idmap[,1]) #set to numeric
  
  nl <- sm$n.latent
  nm <- sm$n.manifest
  
  pdf('VisualParMats.pdf')
  ctsem:::corplotmelt(ctsem:::meltcov(cov2cor(cp$MANIFESTcov)),title =  'Residual Correlations')
  ctsem:::corplotmelt(ctsem:::meltcov(cov2cor(cp$DIFFUSIONcov[latents,latents,drop=FALSE])),title =  'Random latent change Correlations')
  ctsem:::corplotmelt(ctsem:::meltcov(cov2cor(cp$asymDIFFUSION[latents,latents,drop=FALSE])),title =  'Asymptotic latent Correlations')
  dr=cp$DRIFT[latents,latents,drop=FALSE]
  dr[diag(nrow(dr))==1] <- -dr[diag(nrow(dr))==1]
  ctsem:::corplotmelt(ctsem:::meltcov(cov2cor(dr %*% t(dr))),title =  'Std. deterministic change relations')
  dev.off()
  
  ### GENERATING OUTPUT ### 
  
  # generate summary of output in textfile
  sink(file = paste0('summary_',name,'.txt'))
  print(s)
  cat("\n")
  print("Z > 1.96 Correlations:")
  print(s$rawpopcorr[which(abs(s$rawpopcorr[, "z"]) > 1.96), ])
  if (sf$ctstanmodel$n.TIpred > 0) {
    cat("\n")
    print("Z > 1.96 Covariates (TIPs):")
    print(s$tipreds[which(abs(s$tipreds[, "z"]) > 1.96), ])
  }
  cat("\n")
  print("popmeans with 0 NOT within range 2.5% - 97.5%")
  print(
    subset(s$popmeans, s$popmeans$`2.5%` > 0 & s$popmeans$`97.5%` > 0 | s$popmeans$`2.5%` < 0 & s$popmeans$`97.5%` < 0)
  )
  cat("\n")
  print('Par matrices')
  print(cp)
  sink()
  
  #Tables output -> makes word document
  library(flextable)
  library(officer)
  sitems=c('popmeans','popsd','tipreds','rawpopcorr')
  ftlist <- lapply(sitems,
    function(x){
      if(!is.null(s[[x]])){
      fti=data.frame(round(s[[x]] ,2))
      fti=data.frame(Parameter=rownames(fti),fti[,colnames(fti)[!colnames(fti) %in% 'X50.']])
      colnames(fti)[1:5] <- c('Parameter','Mean','SD','2.5%','97.5%')
      return(fti)
      } else return(NULL)
    })
  names(ftlist)=sitems
  
  ftshort <- lapply(c('tipreds','rawpopcorr'), function(x){
    if(!is.null(ftlist[[x]])){
    fts <- ftlist[[x]][apply(ftlist[[x]][,c('2.5%','97.5%')],1,function(x) abs(sum(sign(x)))==2),]
  } else return(NULL)})
  names(ftshort) = c('tipreds','rawpopcorr')
  
  ftlist <- c(ftlist,ftshort)
  ftlist <- ftlist[c(1,2,5,6,3,4)]
  names(ftlist) <- c('Pop. mean parameters','Pop. SD parameters','Significant covariate effects',
    'Significant random effect correlations','All covariate effects','All random effect correlations')
  
  ftlist <- lapply(ftlist,function(x){
    if(!is.null(x)){
    x=autofit(font(
      flextable(as.data.frame(x))
      ,fontname = 'Times New Roman'))
  }
    return(x)})
  save_as_docx(values=ftlist,path='tables.docx')
  
  # #how many subjects have max observed number of time points
  # tldatNArm <- data.frame(ydat)
  # tldatNArm <- tldatNArm[apply(ydat[,sm$manifestNames],1,function(x) any(!is.na(x))),]
  # n<-sapply(unique(tldatNArm$id),function(x) length(tldatNArm$id[tldatNArm$id %in% x]))
  # whichsubfull <- unique(data$id)[n %in% c(max(n):(max(n)-30))]
  # whichsubfull <- match(whichsubfull,sf$setup$idmap[,1]) #set to numeric
  # sink(file = paste0('WhichSubFull_',name,'.txt'))
  # print("WHICH SUBJECTS OBSERVED AT ALL TIME POINTS")
  # cat("\n")
  # print("number of subject with one or more time points")
  # print(length(n))
  # print("Most time points for a subject")
  # print(max(n))
  # print("IDs of subjects with all time points")
  # print(whichsubfull)
  # print("Average Timepoints per ID")
  # print(sum(n)/(length(n)))
  # cat("\n")
  # print("Timepoints per ID")
  # print(n)
  # sink()
  
  #print Latex equation 
  #ctsem:::ctModelLatex(sf,folder = './',filename = paste0('tex_fit'),open = FALSE)
  #ctsem:::ctModelLatex(sm,folder = './',filename = paste0('tex_model'),open = FALSE)
  
  # #TI Correlations matrix plot
  #   sub <- ctsem:::ctStanSubjectPars(sf)
  #   pdf('TI_correlations.pdf')
  #   corm=ctsem:::meltcov(cov2cor(cov(
  #     cbind(sf$standata$tipredsdata,sub[1,,order(dimnames(sub)[[3]])]))))
  #   a=ctsem:::corplotmelt(corm, label = 'Corr.')
  #   print(a+ggplot2::labs(title='Individual difference correlations'))
  #   dev.off()
  #   
  
  
  
  #individual difference / ti pred correlations
  
  #subject par correlation quantiles
  library(plyr)
  library(ggplot2)
  sp=ctsem:::ctStanSubjectPars(sf,pointest=FALSE,cores=cores,nsamples = 500)
  tip <- array(rep(sf$standata$tipredsdata,each=dim(sp)[[1]]),dim=c(dim(sp)[1:2],ncol(sf$standata$tipredsdata)))
  spti=array(c(sp, tip), dim = c(dim(sp)[1], dim(sp)[2], dim(sp)[3]+dim(tip)[3]))
  dn=dimnames(sp)
  dn[[3]] <- c(dn[[3]],sf$ctstanmodelbase$TIpredNames)
  dimnames(spti) <- dn
  
  cornames <- matrix(
    paste0(dimnames(spti)[[3]],'_',rep(dimnames(spti)[[3]],each=length(dimnames(spti)[[3]]))),
    nrow=length(dimnames(spti)[[3]]))
  
  cornames <- cornames[lower.tri(cornames)]
  cm=aaply(spti,1,cor,.drop=FALSE)
  qc=aaply(cm,c(3,2),function(x) round(c(Mean=mean(x),SD=sd(x),
    `2.5%`=quantile(x,probs=c(.025)),`50%`=quantile(x,probs=c(.5)),
    `97.5%`=quantile(x,probs=c(.975)),z=mean(x)/sd(x)),2),.drop=FALSE)
  
  qcsig <- apply(qc[,,c("2.5%.2.5%","97.5%.97.5%"),drop=FALSE],c(1,2),function(x) abs(sum(sign(x)))==2)
  cors <- as.data.table(qc)
  cors <- data.frame(dcast(cors,'V1+V2~V3'))
  colnames(cors) <- c('Par1','Par2','2.5%','50%','97.5%','Mean','SD','z')
  cors <- cors[,c('Par1','Par2','Mean','SD','2.5%','50%','97.5%','z')]
  cors <- cors[paste0(cors$Par1,'_',cors$Par2) %in% cornames,]
  corssig <- cors[apply(cors,1,function(x) sum(sign(as.numeric(x[3:5])))==3),]
  
  #now do similar to point estimate
  spp=ctsem:::ctStanSubjectPars(sf,pointest=T,cores=1)
  tip <- array(rep(sf$standata$tipredsdata,each=dim(spp)[[1]]),dim=c(dim(spp)[1:2],ncol(sf$standata$tipredsdata)))
  spti=array(c(spp, tip), dim = c(dim(spp)[1], dim(spp)[2], dim(spp)[3]+dim(tip)[3]))
  dn=dimnames(spp)
  dn[[3]] <- c(dn[[3]],sf$ctstanmodelbase$TIpredNames)
  dimnames(spti) <- dn
  cm=cor(spti[1,,])
  
  mcor=ctsem:::meltcov(cm)
  msig=data.frame(matrix(as.numeric(qcsig),nrow=nrow(qcsig)))
  dimnames(msig) <- dimnames(cm)
  mcorsig=ctsem:::meltcov(msig)
  mcorsig$sig <- mcorsig$value
  mcorsig$value <- NULL
  mcor <- merge(mcor,mcorsig)
  mcor$`Sig.` <- as.logical(mcor$sig)
  
  pdf('IndDifCorrelations.pdf')
  ggplot(data=(mcor),
    aes_string(x='Var1',y='Var2',fill=('value')))+ #ifelse(groups,NULL,'Group')))+
    geom_tile( width=1,height=1,colour='black')+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
      midpoint = 0, limits = c(-1,1), space = "Lab", 
      name='Corr')  +
    geom_point(mapping=aes(alpha=Sig.),stroke=0,size=.4)+
    theme_minimal()+ theme(axis.text.y=element_text(size=4),
      axis.text.x = element_text(angle = 90,size=4))
  dev.off()
  
  
  
  
  ### CHECKING MODELFIT ###
  
  # Checkfit function 
  pdf('checkfit.pdf')
  by=sm$timeName
  ctsem:::ctCheckFit(fit = sf,data = T,postpred = T,statepred = F,by = by,breaks=2,covplot = T,smooth=T)
  ctsem:::ctCheckFit(fit = sf,data = T,postpred = T,statepred = T,by = by,breaks=4,covplot = F,smooth=T)
  ctsem:::ctCheckFit(fit = sf,data = T,postpred = T,statepred = F,by = by,breaks=4,covplot = T,smooth=T)
  ctsem:::ctCheckFit(fit = sf,data = T,postpred = T,statepred = T,by = by,breaks=2,covplot = T,lag=1,smooth=T)
  dev.off()
  
  
  
  
  if(sf$ctstanmodelbase$n.TIpred > 0){
  #TIP Effect expectations 2
  ktifull=ctsem:::ctKalmanTIP(sf,kalmanvec='yprior',subject = 1,
    #elementNames=c(sf$ctstanmodelbase$latentNames),
    plot=T,polygonsteps=0)
  
  ktifull=ktifull+
    scale_colour_manual(values=c(1,rep(2:((sm$n.TIpred+1)*2),each=2)))
  
  lktifull=ctsem:::ctKalmanTIP(sf,kalmanvec='etaprior',subject = 1,
    plot=T,polygonsteps=0)
  
  lktifull=lktifull+
    scale_colour_manual(values=c(1,rep(2:((1+sm$n.TIpred)*2),each=2)))
  
  pdf('tieffects_expectations.pdf')
  print(lktifull)
  print(ktifull)
  for(ti in sm$TIpredNames){
    ctsem:::ctKalmanTIP(sf,kalmanvec='etaprior',tipreds = ti,subject = whichsubfull[1],
      plot=T,polygonsteps=0)
    ctsem:::ctKalmanTIP(sf,kalmanvec='yprior',tipreds = ti,subject =  whichsubfull[1],
      #elementNames=c('T1asfa','T1ashr','T1itfa','T1ithr','pepdiff'),
      elementNames=c(sf$ctstanmodelbase$manifestNames),
      plot=T,polygonsteps=0)
  }
  dev.off()
  rm(lktifull);rm(ktifull);
  }
  
  
  
  #subject expectation plots 
  k<-ctKalman(sf,subjects = whichsubfull)
  krem<-ctKalman(sf, subjects = whichsubfull,removeObs = TRUE)
  
  pdf('subjectexpectations.pdf')
  plot(krem,polygonsteps=F, kalmanvec='etaprior')
  plot(k,polygonsteps=F, kalmanvec='etasmooth')
  plot(krem,polygonsteps=F, kalmanvec='yprior')
  plot(k,polygonsteps=F, kalmanvec='ysmooth')
  dev.off()
  rm(k);rm(krem)
  
  #Discrete pars plots
  pdf('discretepars.pdf')
  dtpars=ctStanDiscretePars(sf,plot=F,times=times) #,indices=cbind(1:9,i))
  dtparso=ctStanDiscretePars(sf,plot=F,observational = TRUE,times=times) #
  # dtparsc=ctStanDiscretePars(sf,plot=F,observational = TRUE,cov = T,nsamples = 500,times=times) #
  # dtparsoc=ctStanDiscretePars(sf,plot=F,observational = F,cov = T,nsamples = 500,times=times) #
  for(i in 1:sf$ctstanmodelbase$n.latent){
    print(ctStanDiscreteParsPlot(dtpars,indices=cbind(latents,i)))
    print(ctStanDiscreteParsPlot(dtparso,indices=cbind(latents,i)))
    # print(ctStanDiscreteParsPlot(dtparsc,indices=cbind(1:9,i),quantiles = c(.4,.5,.6)))
    # print(ctStanDiscreteParsPlot(dtparsoc,indices=cbind(1:9,i),quantiles = c(.4,.5,.6)))
  }
  dev.off()
 
  
  #paramter posterior plots
  if(sf$standata$nopriors==0){
  pdf('parameter_posterior.pdf')
  ctStanPlotPost(sf,priorwidth = TRUE)
  dev.off()
  }
  
  # #K-fold Cross Validation - OOS entropy score
  #   loo=ctsem:::ctLOO(fit = sf,parallelFolds = F,folds = 10,cores = cores,subjectwise = TRUE,keepfirstobs = F)
  #   save(loo,file='loo.rda')
  #   sink('loo.txt')
  #   print(loo)
  #   sink()  
  
}
}

}
