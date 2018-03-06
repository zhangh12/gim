
mcov <- function(para, para.id, family, data, model, nsample, V, bet0, outcome){
  
  if(family == 'gaussian'){
    mat <- cov.lm(para, para.id, data, model, nsample, V, bet0, outcome)
  }
  
  if(family == 'binomial'){
    mat <- cov.lo(para, para.id, data, model, nsample, V, bet0, outcome)
  }
  
  mat
  
}

