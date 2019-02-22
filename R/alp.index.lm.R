
# return index of alp, model coefficients in the linear model
# error term is ignored
alp.index.lm <- function(map, i){
  
  # no alp except a tau (error term)
  # we always assume that tau is not given
  # so length(map$alp[[i]]) > 0
  if(length(map$alp[[i]]) == 1){
    return(NULL)
  }
  
  map$alp[[i]][-1]
  
}

