
init <- function(formula, family, int, model, nsample){
  
  if(!(family %in% c('gaussian', 'binomial', 'case-control'))){
    stop('family should be \'gaussian\', \'binomial\', or \'case-control\'.')
  }
  
  if(family == 'gaussian'){
    ini <- init.lm(formula, int, model, nsample)
  }
  
  if(family == 'binomial'){
    ini <- init.lo(formula, int, model, nsample)
  }
  
  ini
  
}
