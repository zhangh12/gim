
mcov <- function(para, para.id, family, data, model, nsample, outcome){
  
  if(family == 'gaussian'){
    mat <- cov.lm(para, para.id, data, model, nsample, outcome)
  }
  
  if(family == 'binomial'){
    mat <- cov.lo(para, para.id, data, model, nsample, outcome)
  }
  
  mat
  
}

