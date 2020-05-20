

foo <- function(lam, g, deriv){
  
  n <- nrow(g)
  pr <- 1/as.vector(1 + g %*% lam)
  
  if(deriv == 0){
    return( sum(log(pr))/n )
  }
  
  if(deriv == 1){
    return( -as.vector(t(g) %*% pr)/n )
  }
  
  if(deriv == 2){
    return( t(g) %*% (g * pr^2)/n )
  }
  
}

inner.loop <- function(para, map, ref, sample.info){
  
  ncase <- sample.info$ncase
  nctrl <- sample.info$nctrl
  n <- nrow(ref)
  
  lam0 <- para[map$lam]
  
  tilt <- tilt.cc(para, map, ref)
  Delta <- tilt$Delta
  delta <- tilt$delta
  
  g <- gfunction.cc(para, map, ref, Delta, delta, ncase, nctrl)
  
  niter <- 100
  tol <- 1e-6
  fn <- foo(lam0, g, 0)
  
  for(i in seq(niter)){
    s <- foo(lam0, g, 1)
    if(max(abs(s)) < tol){
      # converged
      break
    }
    
    h <- foo(lam0, g, 2)
    v <- eigen(h)$values
    
    d <- solve(h, -s)
    
    cat(min(v), '\t', max(abs(s)), '\n')
    
    lam1 <- NULL
    for(j in 1:10){
      tau <- .5^(j - 1)
      
      tmp <- as.vector(1 + g %*% (lam0 + tau * d))
      if(all(tmp > 1/n)){
        # convex hull
        lam1 <- lam0 + tau * d
        break
      }else{
        lam1 <- NULL
      }
    }
    
    if(is.null(lam1)){
      stop('no move')
    }
    lam0 <- lam1
    rm(lam1)
    fn <- c(fn, foo(lam0, g, 0))
    plot(fn, pch = 20, main = 'Inner Loop')
  }
  
  mean(1/(1+g%*%lam0))
  s <- foo(lam0, g, 1)
  h <- foo(lam0, g, 2)
  min.ev <- min(eigen(h)$value)
  if(min.ev < 0){
    warning(paste0('Inner loop may stop at a saddle point: min ev = ', min.ev))
  }
  
  if(max(abs(s)) > tol){
    warning(paste0('Inner loop does not stop at a stationary points: max grad = ', max(abs(s))))
  }
  
  #message(paste0('Inner loop converges with min ev = ', min.ev, '; max grad = ', max(abs(s))))
  
  lam0
  
}
