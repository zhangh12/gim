

gim.default <- function(formula, int, model, nsample, outcome = 'y'){
  
  ini <- init.lm(formula, int, model, nsample)
  para <- ini$para
  para.id <- ini$para.id
  bet0 <- ini$bet0
  
  V <- Sigma0.lm(para, para.id, int, model, nsample)
  fit <- NR.lm(para, para.id, int, V, bet0)
  fit$vcov <- cov.lm(fit$coefficients, para.id, int, model, nsample)
  
  fit <- reorganize(fit, para.id)
  
  fit$call <- match.call()
  class(fit) <- 'gim'
  
  fit
  
}
