
mcov <- function(para, map, family, data, ref, model, nsample, V, bet0, outcome){
  
  if(family == 'gaussian'){
    mat <- cov.lm(para, map, data, ref, model, nsample, V, bet0, outcome)
  }
  
  if(family == 'binomial'){
    mat <- cov.lo(para, map, data, ref, model, nsample, V, bet0, outcome)
  }
  
  mat
  
}

