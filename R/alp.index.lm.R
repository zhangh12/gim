
# return index of alp, model coefficients in the linear model
alp.index.lm <- function(map, i){
  
  # a logistic regression model, all coef given (inc intercept), no alp for this model
  if(any(is.na(map$alp[[i]]))){
    return(NULL)
  }
  
  map$alp[[i]]
  
}

