

optimal.Sigma0 <- function(para, para.id, family, int, model, nsample, outcome){
  
  if(family == 'gaussian'){
    V <- Sigma0.lm(para, para.id, int, model, nsample, outcome)
  }
  
  if(family == 'binomial'){
    V <- Sigma0.lo(para, para.id, int, model, nsample, outcome)
  }
  
  V
  
}

