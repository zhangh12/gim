

convex.hull <- function(para, map, family, ref, sample.info, step, p){
  
  if(family == 'gaussian'){
    step <- check.convex.hull.lm(para, map, ref, step, p)
  }
  
  if(family == 'binomial'){
    
  }
  
  if(family == 'case-control'){
    step <- check.convex.hull.cc(para, map, ref, sample.info, step, p)
  }
  
  if(family == 'cml'){
    
  }
  
  return(step)
  
}
