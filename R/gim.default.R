

gim.default <- function(formula, family, int, model, nsample, outcome = 'y'){
  
  ini <- init(formula, family, int, model, nsample)
  para <- ini$para
  para.id <- ini$para.id
  bet0 <- ini$bet0
  
  V <- optimal.Sigma0(para, para.id, family, int, model, nsample, outcome)
  fit <- NR(para, para.id, family, int, V, bet0, outcome)
  fit$vcov <- mcov(fit$coefficients, para.id, family, int, model, nsample, outcome)
  
  fit <- reorganize(fit, para.id, family)
  
  fit$call <- match.call()
  class(fit) <- 'gim'
  
  fit
  
}
