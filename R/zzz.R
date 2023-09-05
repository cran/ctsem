# # # .onLoad <- function(libname, pkgname) {
#   checkm<-ctModel(
#     type='stanct',
#     n.latent=2,n.TDpred=1,n.TIpred=1,n.manifest=2,
#     MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
#     MANIFESTMEANS=0,
#     DRIFT=c('dr1','dr12','dr21||||TI1','dr22'),
#     DIFFUSION=c('diff11',0,'diff21','diff22||||TI1'),
#     CINT=matrix(c('cint1||||TI1','cint2||||TI1'),ncol=1),
#     LAMBDA=diag(2),tipredDefault=FALSE)
# 
#   ctstantestfit<-ctStanFit(ctsem::ctstantestdat,checkm,cores=1,inits=0,
#     optimize = TRUE,optimcontrol=list(finishsamples=20,stochastic=F,carefulfit=F),priors=TRUE)
# 
#   ctstantestfit <- ctStanGenerateFromFit(ctstantestfit,nsamples = 20,fullposterior = TRUE,cores=1)
#   ("ctstantestfit", ctstantestfit, envir = asNamespace('ctsem'))
# # # }
