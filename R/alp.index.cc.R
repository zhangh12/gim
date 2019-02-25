
# return index of alp, log OR of all variables, including intercept (if any)
alp.index.cc <- function(map, i){
  
  # a logistic regression model, all coef given (inc intercept), no alp for this model
  if(any(is.na(map$alp[[i]]))){
    return(NULL)
  }
  
  map$alp[[i]]
  
}

