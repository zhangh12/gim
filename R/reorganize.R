

reorganize <- function(fit, map, family){
  
  if(family == 'gaussian'){
    fit <- reorganize.lm(fit, map)
  }
  
  if(family == 'binomial'){
    fit <- reorganize.lo(fit, map)
  }
  
  fit
  
}
