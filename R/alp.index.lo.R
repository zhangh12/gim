
alp.index.lo <- function(id.alp, i){
  
  st <- id.alp$start[i]
  ed <- id.alp$end[i]
  
  if(any(is.na(c(st, ed)))){ # a logistic regression model, all coef given, no alp for this model
    return(NULL)
  }
  
  if(st > ed){
    stop('debug alp.index.lo')
  }
  
  st:ed
  
}

