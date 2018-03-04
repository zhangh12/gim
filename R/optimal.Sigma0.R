

optimal.Sigma0 <- function(para, para.id, family, data, model, nsample, outcome){
  
  if(family == 'gaussian'){
    V <- Sigma0.lm(para, para.id, data, model, nsample, outcome)
  }
  
  if(family == 'binomial'){
    V <- Sigma0.lo(para, para.id, data, model, nsample, outcome)
  }
  
  V
  
}

