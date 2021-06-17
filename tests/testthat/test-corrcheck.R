if(identical(Sys.getenv("NOT_CRAN"), "true")& .Machine$sizeof.pointer != 4){
  library(ctsem)
  library(testthat)
  cores=2
  
  # context("corrCheck")
  
  #anomauth
  test_that("corrCheck", {
    set.seed(1)
    gm <- ctModel(LAMBDA = diag(1,3),DRIFT=diag(-1,3),
      T0VAR=matrix(c(5,-5,-5,0,1,-1,0,0,2),3,3),
      MANIFESTTRAITVAR=matrix(c(2,-1,-1, 0,1,1,0,0,2),3,3),
      DIFFUSION=matrix(c(2,1,1,0,4,-2,0,0,2),3,3),Tpoints=30)
    d <- ctGenerate(ctmodelobj = gm,n.subjects = 600,burnin = 0)
    
    m <- ctModel(LAMBDA = diag(1,3),DRIFT=diag(-1,3),type='stanct',
      # MANIFESTMEANS = 0,
      MANIFESTVAR = 0)
    
    f <- ctStanFit(datalong = d,ctstanmodel = m,nopriors=F)
    
    p <- ctStanContinuousPars(f)
    
    diffcov <- p$DIFFUSIONcov
    diffcor <- cov2cor(p$DIFFUSIONcov)
    ediffcov <- tcrossprod(gm$DIFFUSION)
    ediffcor <- cov2cor(tcrossprod(gm$DIFFUSION))
    
    s=summary(f)
    s$rawpopcorr
    f$stanfit$transformedparsfull$rawpopcorr[1,,]
    tcrossprod(gm$MANIFESTTRAITVAR)
    cov2cor(tcrossprod(gm$MANIFESTTRAITVAR))
    
    #check diagonal of 1's for corr
    testthat::expect_equivalent(diag(f$stanfit$transformedparsfull$rawpopcorr[1,,]),
      rep(1,nrow(f$stanfit$transformedparsfull$rawpopcorr[1,,])),tol=1e-5)
    
    #cov check
    testthat::expect_equivalent(f$stanfit$transformedparsfull$popcov[1,4:6,4:6],
      tcrossprod(gm$MANIFESTTRAITVAR),tol=.5)
    
    #cor check
    testthat::expect_equivalent(f$stanfit$transformedparsfull$rawpopcorr[1,4:6,4:6],
      cov2cor(tcrossprod(gm$MANIFESTTRAITVAR)),tol=1e-1)
    
    
  })
}

