
compute.hess <- function(para, map, family, data, ref, inv.V, bet0, sample.info, outcome, type){
  
  if(family == 'gaussian'){
    h0 <- hess.lm(para, map, data, ref, inv.V, bet0, outcome)
  }
  
  if(family == 'binomial'){
    h0 <- hess.lo(para, map, data, ref, inv.V, bet0, outcome)
  }
  
  if(family == 'case-control'){
    if(type == 'cc-ref'){
      h0 <- hess.ccr(para, map, data, ref, inv.V, bet0, sample.info, outcome)
    }else{
      h0 <- hess.cc(para, map, data, ref, inv.V, bet0, sample.info, outcome)
    }
  }
  
  if(family == 'cml'){
    h0 <- hess.cml(para, map, data, ref, inv.V, bet0, sample.info, outcome)
  }
  
  h0
  
}
