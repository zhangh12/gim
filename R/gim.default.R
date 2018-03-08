

gim.default <- function(formula, family, data, model, nsample){
  
  nsample <- as.matrix(nsample)
  formula <- as.formula(formula)
  
  fp <- formula.parse(formula, model, data)
  model <- fp$model
  data <- fp$data
  outcome <- fp$outcome
  
  ini <- init(formula, family, data, model, nsample)
  para <- ini$para
  para.id <- ini$para.id
  bet0 <- ini$bet0
  
  message('Running Newton-Raphson algorithm on first stage...')
  V <- optimal.Sigma0(para, para.id, family, data, model, nsample, outcome)
  para <- NR(para, para.id, family, data, V, bet0, outcome)$coefficients
  
  message('Running Newton-Raphson algorithm on second stage...')
  V <- optimal.Sigma0(para, para.id, family, data, model, nsample, outcome)
  fit <- NR(para, para.id, family, data, V, bet0, outcome)
  
  fit$vcov <- mcov(fit$coefficients, para.id, family, data, model, nsample, V, bet0, outcome)
  
  fit <- reorganize(fit, para.id, family)
  
  fit$call <- match.call()
  class(fit) <- 'gim'
  
  fit
  
}
