
check.convex.hull.lm <- function(para, map, ref, step, p){
  
  ret <- NULL
  for(a in step){
    lam <- (para + a * p)[map$lam]
    g <- gfunction.lm(para + a * p, map, ref)
    
    v <- as.vector(1+g %*% lam)
    
    n <- nrow(ref)
    if(all(v > 1/n)){
      ret <- c(ret, a)
    }
  }
  
  if(is.null(ret)){
    stop('convex hull fails')
  }
  
  return(ret)
  
}
