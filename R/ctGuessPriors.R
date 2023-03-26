

ctGuessPriors <- function(dat, ctm){
  message('Roughly estimating workable priors, assuming relatively stationary data, factor loadings ~ 1, and simple model structure.')
  l=listOfMatrices(ctm$pars)
  
  lambda <- l$LAMBDA
  lambda[is.na(as.numeric(lambda))] <- 1
  ilambda <- MASS::ginv(lambda)
  manifests <- dat[,ctm$manifestNames,drop=FALSE]
  latents <- matrix(apply(manifests,1,function(x) matrix(t(ilambda) %*% x,nrow=ncol(lambda))),ncol=ncol(lambda))
  t0latents <-matrix(latents[match(unique(dat[,ctm$subjectIDname]),dat[,ctm$subjectIDname]),,drop=FALSE],ncol=ncol(latents))
  
  t0priors <-  paste0(apply(t0latents,2,function(x) round(mean(x,na.rm=TRUE),3)), #mean of first observed data
    ' + param * ', apply(t0latents,2,function(x) round(sd(x,na.rm=TRUE),3)))
  
  manifestsd <- apply(manifests,2,function(x) round(sd(x,na.rm=TRUE),3))
  
  mvarpriors <-  paste0(manifestsd,' * log1p(exp(param*2-2))')
  
  manifestmeanpriors <-  paste0(apply(manifests,2,function(x) round(mean(x,na.rm=TRUE),3)), #mean of first observed data
    ' + param * ', manifestsd *2 )
  
  if(any(is.na(ctm$pars$value[ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$row <= ncol(latents)]))){ #t0means
    ctm$pars$transform[ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$row <= ncol(latents) & is.na(ctm$pars$value)] <- 
      t0priors[ctm$pars$row[ctm$pars$matrix %in% 'T0MEANS' & ctm$pars$row <= ncol(latents)& is.na(ctm$pars$value)]]
  }
  
  if(any(is.na(ctm$pars$value[ctm$pars$matrix %in% 'MANIFESTMEANS']))){ 
    ctm$pars$transform[ctm$pars$matrix %in% 'MANIFESTMEANS' & is.na(ctm$pars$value)] <- 
      manifestmeanpriors[ctm$pars$row[ctm$pars$matrix %in% 'MANIFESTMEANS' & is.na(ctm$pars$value)]]
  }
  
  if(any(is.na(ctm$pars$value[ctm$pars$matrix %in% 'MANIFESTVAR']))){ 
    ctm$pars$transform[ctm$pars$matrix %in% 'MANIFESTVAR' & is.na(ctm$pars$value)] <- 
      mvarpriors[ctm$pars$row[ctm$pars$matrix %in% 'MANIFESTVAR' & is.na(ctm$pars$value)]]
  }
  
  return(ctm)
}


