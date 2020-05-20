

compute.score <- function(para, map, family, data, ref, inv.V, bet0, sample.info, outcome, type){
  
  if(family == 'gaussian'){
    s0 <- score.lm(para, map, data, ref, inv.V, bet0, outcome)
    #s1 <- grad(obj.lm, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome)
  }
  
  if(family == 'binomial'){
    s0 <- score.lo(para, map, data, ref, inv.V, bet0, outcome)
    #s1 <- grad(obj.lo, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome)
  }
  
  if(family == 'case-control'){
    if(type == 'cc-ref'){
      s0 <- score.ccr(para, map, data, ref, inv.V, bet0, sample.info, outcome)
      #s1 <- grad(obj.ccr, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, sample.info = sample.info, outcome = outcome)
    }else{
      s0 <- score.cc(para, map, data, ref, inv.V, bet0, sample.info, outcome)
      #s1 <- grad(obj.cc, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, sample.info = sample.info, outcome = outcome)
    }
    
  }
  
  if(family == 'cml'){
    s0 <- score.cml(para, map, data, ref, inv.V, bet0, sample.info, outcome)
    #s1 <- grad(obj.cml, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, sample.info = sample.info, outcome = outcome)
  }
  
  s0
  
}
