
# return covariance of auxiliary estimates, e.g. V/N (NOT V) in the paper
optimal.Sigma0 <- function(para, para.id, family, ref, model, nsample, outcome){
  
  if(family == 'gaussian'){
    V <- Sigma0.lm(para, para.id, ref, model, nsample, outcome)
  }
  
  if(family == 'binomial'){
    V <- Sigma0.lo(para, para.id, ref, model, nsample, outcome)
  }
  
  V
  
}

