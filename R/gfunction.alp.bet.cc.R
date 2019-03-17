

gfunction.alp.bet.cc <- function(para, map, ref, Delta, delta, ncase, nctrl, xi, pr){
  
  nmodel <- length(map$bet)
  
  g.alp.bet <- list()
  
  n <- nrow(ref)
  
  const <- list()
  for(i in 1:nmodel){
    id <- c(alp.index.cc(map, i), map$bet[[i]])
    gam <- para[id]
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    
    rho.i <- ncase[i, i] / nctrl[i, i]
    const[[i]] <- -rx * (delta[, i] * (1 - rho.i * delta[, i]) * (1 + rho.i * Delta) / (1 + rho.i * delta[, i])^3)
  }
  
  nlam <- max(map$lam)
  offset <- max(map$the)
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  for(i in 1:nmodel){
    id.a <- alp.index.cc(map, i)
    if(is.null(id.a)){
      next
    }
    
    id.b <- map$bet[[i]]
    
    id <- c(id.a, id.b)
    for(j in id.a){
      fxj <- ref[, names(para)[j]]
      tmp <- const[[i]] * fxj
      for(l in id.b){
        fxl <- ref[, names(para)[l]]
        gt <- matrix(0, nrow = n, ncol = nlam - 1)
        gt[, id - offset] <- tmp * fxl
        #g.alp.bet[[foo(j,l)]] <- gt
        g.alp.bet[[foo(j,l)]] <- t(gt %*% xi) %*% pr
      }
    }
  }
  
  if(length(g.alp.bet) == 0){
    g.alp.bet <- NULL
  }
  
  g.alp.bet
  
}


