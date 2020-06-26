

gfunction.the2.lo <- function(para, map, ref, Delta, lam, pr){
  
  nmodel <- length(map$bet)
  
  the <- para[map$the]
  
  g.the2 <- list()
  
  nthe <- length(the)
  n <- nrow(ref)
  
  const <- list()
  tmp <- Delta * (1 - Delta) / (1 + Delta)^3
  for(i in 1:nmodel){
    id <- c(alp.index.lo(map, i), map$bet[[i]])
    gam <- para[id]
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    
    const[[i]] <- rx * tmp
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
      gt <- matrix(0, nrow = n, ncol = nlam)
      for(i in 1:nmodel){
        id <- c(alp.index.lo(map, i), map$bet[[i]])
        gt[, id - offset] <- const[[i]] * fxj * fxl
      }
      
      g.the2[[foo(j,l)]] <- t(gt %*% lam) %*% pr
      rm(gt)
    }
  }
  
  g.the2
  
}


