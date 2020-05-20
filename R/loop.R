
loop <- function(para, map, family, data, ref, V, bet0, sample.info, outcome, type, silent){
  
  inv.V <- solve(V)
  np <- length(para)
  n <- nrow(ref)
  
  x0 <- para
  fn <- compute.obj(x0, map, family, data, ref, inv.V, bet0, sample.info, outcome, type)
  
  niter <- 50
  tol <- 1e-6
  
  for(i in seq(niter)){
    x0[map$lam] <- inner.loop(x0, map, ref, sample.info) ## lambda(mu)
    s <- compute.score(x0, map, family, data, ref, inv.V, bet0, sample.info, outcome, type)
    h <- compute.hess(x0, map, family, data, ref, inv.V, bet0, sample.info, outcome, type)
    
    s <- s[-map$lam]
    h <- h[-map$lam, -map$lam] - h[-map$lam, map$lam] %*% solve(h[map$lam, map$lam]) %*% h[map$lam, -map$lam]
    
    ei <- eigen(h)
    if(any(ei$values >= -1e-2)){
      message('modifying indefinite hess')
      ev <- ei$values
      ev <- ev * (ev < -1e-2) - 0.01 * (ev >= -1e-2)
      h <- ei$vectors %*% diag(ev) %*% t(ei$vectors)
    }
    
    max.ev <- max(eigen(h)$values)
    deriv1 <- max(abs(s))
    max.fn <- tail(fn, 1)
    
    if(deriv1 < tol && max.ev < 0){
      ## converge
      break
    }
    
    d <- solve(h, s)
    x1 <- NULL
    for(j in 1:20){
      tau <- .5^(j - 1)
      x1 <- x0
      x1[-map$lam] <- x1[-map$lam] - tau * d
      x1[map$lam] <- inner.loop(x1, map, ref, sample.info)
      
      ncase <- sample.info$ncase
      nctrl <- sample.info$nctrl
      tilt <- tilt.cc(x1, map, ref)
      Delta <- tilt$Delta
      delta <- tilt$delta
      
      g <- gfunction.cc(x1, map, ref, Delta, delta, ncase, nctrl)
      tmp <- as.vector(1 + g %*% x1[map$lam])
      fn1 <- compute.obj(x1, map, family, data, ref, inv.V, bet0, sample.info, outcome, type)
      if(all(tmp > 1/n) && fn1 > max.fn){
        # a valid move in outer loop
        break
      }else{
        x1 <- NULL
      }
    }
    
    if(is.null(x1)){
      stop('cannot find search direction in outer loop')
    }
    
    x0 <- x1
    fn <- c(fn, fn1)
    plot(fn, pch = 20, main = 'Outer Loop')
  }
  
  deriv1 <- formatC(deriv1, digits = 1, format = 'e')
  max.ev <- formatC(max.ev, digits = 1, format = 'e')
  max.fn <- formatC(tail(fn, 1), digits = 1, format = 'e')
  
  if(max(abs(s)) < tol && max.ev < 0){
    message(paste0('outer loop converges. grad = ', deriv1, ', max.ev = ', max.ev, ', max.fn = ', max.fn))
  }else{
    if(max(abs(s)) >= tol){
      message(paste0('outer loop cannot find a stationary point. grad = ', deriv1, ', max.ev = ', max.ev))
    }else{
      message(paste0('outer loop stops at a saddle point. grad = ', deriv1, ', max.ev = ', max.ev))
    }
  }
  
  x0
  
}
