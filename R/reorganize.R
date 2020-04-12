

reorganize <- function(fit, map, family, type){
  
  if(family == 'gaussian'){
    fit <- reorganize.lm(fit, map)
  }
  
  if(family == 'binomial'){
    fit <- reorganize.lo(fit, map)
  }
  
  if(family == 'case-control'){
    if(type == 'cc-ref'){
      fit <- reorganize.ccr(fit, map)
    }else{
      fit <- reorganize.cc(fit, map)
    }
  }
  
  if(family == 'cml'){
    fit <- reorganize.cml(fit, map)
  }
  
  fit
  
}
