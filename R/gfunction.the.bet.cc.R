

gfunction.the.bet.cc <- function(para, map, ref, Delta, delta, ncase, nctrl, xi, pr){
  
  nmodel <- length(map$bet)
  
  the <- para[map$the]
  
  g.the.bet <- list()
  
  nthe <- length(the)
  n <- nrow(ref)
  
  const <- list()
  for(i in 1:nmodel){
    id <- c(alp.index.cc(map, i), map$bet[[i]])
    gam <- para[id]
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    
    rho.i <- ncase[i, i] / nctrl[i, i]
    const[[i]] <- -rx * (Delta * rho.i * delta[, i] / (1 + rho.i * delta[, i])^2)
  }
  
  nlam <- max(map$lam)
  offset <- max(map$the)
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  for(j in 1:nthe){
    fxj <- ref[, names(the)[j]]
    for(i in 1:nmodel){
      id.b <- map$bet[[i]]
      id <- c(alp.index.cc(map, i), map$bet[[i]])
      tmp <- const[[i]] * fxj
      for(l in id.b){
        fxl <- ref[, names(para)[l]]
        gt <- matrix(0, nrow = n, ncol = nlam - 1)
        gt[, id - offset] <- tmp * fxl
        #g.the.bet[[foo(j,l)]] <- gt
        g.the.bet[[foo(j,l)]] <- t(gt %*% xi) %*% pr
        rm(gt)
      }
    }
  }
  
  g.the.bet
  
}


