

reorganize <- function(fit, para.id, family){
  
  if(family == 'gaussian'){
    fit <- reorganize.lm(fit, para.id)
  }
  
  if(family == 'binomial'){
    fit <- reorganize.lo(fit, para.id)
  }
  
  fit
  
}
