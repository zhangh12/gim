
mcov <- function(para, map, family, data, ref, model, sample.info, V, bet0, outcome){
  
  if(family == 'gaussian'){
    mat <- cov.lm(para, map, data, ref, model, sample.info, V, bet0, outcome)
  }
  
  if(family == 'binomial'){
    mat <- cov.lo(para, map, data, ref, model, sample.info, V, bet0, outcome)
  }
  
  if(family == 'case-control'){
    mat <- cov.cc(para, map, data, ref, sample.info, V, bet0, outcome)
  }
  
  mat
  
}

