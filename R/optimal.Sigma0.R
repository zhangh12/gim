
# return covariance of auxiliary estimates, e.g. V/N (NOT V) in the paper
optimal.Sigma0 <- function(para, map, family, ref, model, sample.info, pr0, Delta, outcome, type, bet0 = NULL){
  
  if(family == 'gaussian'){
    V <- Sigma0.lm(para, map, ref, model, sample.info, outcome)
  }
  
  if(family == 'binomial'){
    V <- Sigma0.lo(para, map, ref, model, sample.info, outcome)
  }
  
  if(family == 'case-control' || family == 'cml'){
    if(family == 'cml'){
      para[map$all.bet] <- bet0
      names(para)[map$all.bet] <- names(bet0)
    }
    if(type == 'cc-ref'){
      V <- Sigma0.ccr(para, map, ref, model, sample.info, pr0, Delta, outcome)
    }else{
      V <- Sigma0.cc(para, map, ref, model, sample.info, pr0, Delta, outcome)
    }
  }
  
  V
  
}

