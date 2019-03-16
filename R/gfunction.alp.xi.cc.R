

gfunction.alp.xi.cc <- function(g.alp, xi){
  
  if(is.null(g.alp)){
    return(NULL)
  }
  
  n <- nrow(g.alp[[1]])
  nalp <- length(g.alp)
  
  g.alp.xi <- matrix(NA, nrow = n, ncol = nalp)
  k <- 0
  for(m in g.alp){
    k <- k + 1
    g.alp.xi[, k] <- m %*% xi
  }
  
  g.alp.xi
  
}
