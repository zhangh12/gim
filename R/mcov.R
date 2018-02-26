
mcov <- function(para, para.id, family, int, model, nsample, outcome){
  
  if(family == 'gaussian'){
    mat <- cov.lm(para, para.id, int, model, nsample, outcome)
  }
  
  if(family == 'binomial'){
    mat <- cov.lo(para, para.id, int, model, nsample, outcome)
  }
  
  mat
  
}

