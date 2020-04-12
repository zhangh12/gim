
init <- function(fit0, family, data, ref, model, nsample, ncase, nctrl, outcome, type){
  
  if(!(family %in% c('gaussian', 'binomial', 'case-control', 'cml'))){
    stop('family should be \'gaussian\', \'binomial\', or \'case-control\'.')
  }
  
  if(family == 'gaussian'){
    ini <- init.lm(fit0, data, model, nsample)
  }
  
  if(family == 'binomial'){
    ini <- init.lo(fit0, data, model, nsample)
  }
  
  if(family == 'case-control'){
    if(type == 'cc-ref'){
      ini <- init.ccr(fit0, data, ref, model, ncase, nctrl, outcome)
    }else{
      ini <- init.cc(fit0, data, model, ncase, nctrl, outcome)
    }
  }
  
  if(family == 'cml'){
    ini <- init.cml(fit0, data, model, ncase, nctrl, outcome)
  }
  
  ini
  
}
