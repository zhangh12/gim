

gfunction.the.cc <- function(para, map, ref, Delta, delta, ncase, nctrl){
  
  nmodel <- length(map$bet)
  
  the <- para[map$the]
  
  g.the <- list()
  
  nthe <- length(the)
  n <- nrow(ref)
  
  const <- list()
  for(i in 1:nmodel){
    id <- c(alp.index.cc(map, i), map$bet[[i]])
    gam <- para[id]
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    
    rho.i <- ncase[i, i] / nctrl[i, i]
    const[[i]] <- rx * (Delta / (1 + rho.i * delta[, i]))
  }
  
  nlam <- max(map$lam)
  #offset <- max(map$the) - 1
  offset <- max(map$the)
  for(j in 1:nthe){
    gt <- matrix(0, nrow = n, ncol = nlam - 1)
    fx0 <- ref[, names(the)[j]]
    #gt[, 1] <- fx0 * Delta
    for(i in 1:nmodel){
      id <- c(alp.index.cc(map, i), map$bet[[i]])
      gt[, id - offset] <- const[[i]] * fx0
    }
    g.the[[j]] <- gt
    rm(gt)
  }
  
  g.the
  
}


