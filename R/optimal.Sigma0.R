
# return covariance of auxiliary estimates, e.g. V/N (NOT V) in the paper
optimal.Sigma0 <- function(para, map, family, ref, model, sample.info, pr0, Delta, outcome){
  
  if(family == 'gaussian'){
    V <- Sigma0.lm(para, map, ref, model, sample.info, outcome)
  }
  
  if(family == 'binomial'){
    V <- Sigma0.lo(para, map, ref, model, sample.info, outcome)
  }
  
  if(family == 'case-control'){
    V <- Sigma0.cc(para, map, ref, model, sample.info, pr0, Delta, outcome)
  }
  
  V
  
}

