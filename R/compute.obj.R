

compute.obj <- function(para, map, family, data, ref, inv.V, bet0, sample.info, outcome, type){
  
  if(family == 'gaussian'){
    o <- obj.lm(para, map, data, ref, inv.V, bet0, outcome)
  }
  
  if(family == 'binomial'){
    o <- obj.lo(para, map, data, ref, inv.V, bet0, outcome)
  }
  
  if(family == 'case-control'){
    o <- obj.cc(para, map, data, ref, inv.V, bet0, sample.info, outcome)
  }
  
  if(family == 'cml'){
    o <- obj.cml(para, map, data, ref, inv.V, bet0, sample.info, outcome)
  }
  
  o
  
}
