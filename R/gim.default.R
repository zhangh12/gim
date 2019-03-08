

gim.default <- function(formula, family, data, model, 
                        nsample = NULL, ncase = NULL, nctrl = NULL, 
                        ref = NULL, niter = 2){
  
  data <- clean.data(data)
  
  fp <- formula.parse(formula, model, data, ref)
  model <- fp$model
  data <- fp$data
  ref <- fp$ref
  outcome <- fp$outcome
  formula <- fp$formula
  
  ini <- init(formula, family, data, model, nsample, ncase, nctrl, outcome)
  para <- ini$para
  map <- ini$map
  bet0 <- ini$bet0
  pr0 <- ini$pr0
  Delta <- ini$Delta
  sample.info <- ini$sample.info
  tau <- ini$tau
  
  if(niter < 2){
    stop('niter should at least be 2')
  }
  
  while(niter > 0){
    #message('Running Newton-Raphson algorithm on first stage...')
    V <- optimal.Sigma0(para, map, family, ref, model, sample.info, pr0, Delta, tau, outcome)
    fit <- NR(para, map, family, data, ref, V, bet0, sample.info, outcome)
    para <- fit$coefficients
    niter <- niter - 1
  }
  
  fit$vcov <- mcov(fit$coefficients, map, family, data, ref, model, sample.info, V, bet0, outcome)
  
  fit <- reorganize(fit, map, family)
  
  fit$call <- match.call()
  fit$V.bet <- V
  class(fit) <- 'gim'
  
  fit
  
}
