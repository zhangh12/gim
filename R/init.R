
init <- function(formula, family, data, model, nsample, ncase, nctrl, outcome){
  
  if(!(family %in% c('gaussian', 'binomial', 'case-control'))){
    stop('family should be \'gaussian\', \'binomial\', or \'case-control\'.')
  }
  
  if(family == 'gaussian'){
    ini <- init.lm(formula, data, model, nsample)
  }
  
  if(family == 'binomial'){
    ini <- init.lo(formula, data, model, nsample)
  }
  
  if(family == 'case-control'){
    ini <- init.cc(formula, data, model, ncase, nctrl, outcome)
  }
  
  ini
  
}
