

gim.default <- function(formula, family, data, model, 
                        nsample = NULL, ncase = NULL, nctrl = NULL, 
                        ref = NULL, ...){
  
  argu <- list(...)
  if(is.null(argu$niter)){
    niter <- 100
  }else{
    niter <- argu$niter
    if(niter < 2){
      stop('niter should at least be 2')
    }
  }
  
  if(is.null(argu$tol)){
    tol <- 1e-4
  }else{
    tol <- argu$tol
  }
  
  if(is.null(argu$silent)){
    silent <- TRUE
  }else{
    silent <- FALSE
  }
  
  data <- clean.data(data)
  
  fp <- formula.parse(formula, family, data, model, ref)
  model <- fp$model
  data <- fp$data
  ref <- fp$ref
  outcome <- fp$outcome
  fit0 <- fp$fit0
  type <- fp$type
  
  ini <- init(fit0, family, data, ref, model, nsample, ncase, nctrl, outcome, type)
  para <- ini$para
  map <- ini$map
  bet0 <- ini$bet0
  pr0 <- ini$pr0
  Delta <- ini$Delta
  sample.info <- ini$sample.info
  
  eps <- 1.
  
  while(niter > 0){
    #message('Running Newton-Raphson algorithm on first stage...')
    V <- optimal.Sigma0(para, map, family, ref, model, sample.info, pr0, Delta, outcome, type, bet0)
    if(eps < tol){
      break
    }
    fit <- NR(para, map, family, data, ref, V, bet0, sample.info, outcome, type, silent)
    eps <- max(abs(para - fit$coefficients))
    para <- fit$coefficients
    niter <- niter - 1
  }
  
  fit$vcov <- mcov(para, map, family, data, ref, model, sample.info, V, bet0, outcome, type)
  
  fit <- reorganize(fit, map, family, type)
  
  fit$call <- match.call()
  fit$V.bet <- V
  class(fit) <- 'gim'
  
  fit
  
}
