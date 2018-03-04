

gim.default <- function(formula, family, data, model, nsample, outcome = 'y'){
  
  nsample <- as.matrix(nsample)
  
  ini <- init(formula, family, data, model, nsample)
  para <- ini$para
  para.id <- ini$para.id
  bet0 <- ini$bet0
  
  V <- optimal.Sigma0(para, para.id, family, data, model, nsample, outcome)
  fit <- NR(para, para.id, family, data, V, bet0, outcome)
  fit$vcov <- mcov(fit$coefficients, para.id, family, data, model, nsample, outcome)
  
  fit <- reorganize(fit, para.id, family)
  
  fit$call <- match.call()
  class(fit) <- 'gim'
  
  fit
  
}
