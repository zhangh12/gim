

gfunction.the2.cc <- function(para, map, ref, Delta, delta, ncase, nctrl, xi, pr){
  
  nmodel <- length(map$bet)
  
  the <- para[map$the]
  
  g.the2 <- list()
  
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
  offset <- max(map$the)
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  for(j in 1:nthe){
    fxj <- ref[, names(the)[j]]
    for(l in j:nthe){
      fxl <- ref[, names(the)[l]]
      gt <- matrix(0, nrow = n, ncol = nlam - 1)
      gt[, 1] <- fxj^2 * Delta
      for(i in 1:nmodel){
        id <- c(alp.index.cc(map, i), map$bet[[i]])
        gt[, id - offset] <- const[[i]] * fxj * fxl
      }
      
      # g.the2[[foo(j,l)]] <- gt
      g.the2[[foo(j,l)]] <- t(gt %*% xi) %*% pr
      rm(gt)
    }
  }
  
  g.the2
  
}


