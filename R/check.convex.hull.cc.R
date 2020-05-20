
check.convex.hull.cc <- function(para, map, ref, sample.info, step, p){
  
  ncase <- sample.info$ncase
  nctrl <- sample.info$nctrl
  
  ret <- NULL
  for(a in step){
    lam <- (para + a * p)[map$lam]
    
    tilt <- tilt.cc(para + a * p, map, ref)
    Delta <- tilt$Delta
    delta <- tilt$delta
    
    g <- gfunction.cc(para + a * p, map, ref, Delta, delta, ncase, nctrl)
    
    v <- as.vector(1+g %*% lam)
    n <- nrow(ref)
    print(a)
    if(all(v > 1/n)){
      ret <- c(ret, a)
    }
  }
  
  if(is.null(ret)){
    stop('convex hull fails')
  }
  
  return(ret)
  
}
