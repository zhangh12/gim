
alp.index.lm <- function(id.alp, i){
  
  st <- id.alp$start[i]
  ed <- id.alp$end[i]
  
  if(st > ed){
    stop('debug alp.index.lm')
  }
  
  if(st == ed){ # all coefficients are known for this model, no alp except a tau
    return(NULL)
  }
  
  (st + 1):ed
  
}

