

gfunction.alp2.lo <- function(para, map, ref, delta, lam, pr){
  
  nmodel <- length(map$bet)
  
  g.alp2 <- list()
  
  n <- nrow(ref)
  
  const <- list()
  for(i in 1:nmodel){
    id <- c(alp.index.lo(map, i), map$bet[[i]])
    gam <- para[id]
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    
    const[[i]] <- -rx * (delta[, i] * (1 - delta[, i]) / (1 + delta[, i])^3)
  }
  
  nlam <- max(map$lam)
  offset <- max(map$the)
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  for(i in 1:nmodel){
    id.a <- alp.index.lo(map, i)
    if(is.null(id.a)){
      next
    }
    
    id <- c(id.a, map$bet[[i]])
    for(j in id.a){
      fxj <- ref[, names(para)[j]]
      tmp <- const[[i]] * fxj
      for(l in id.a){
        if(l < j){
          next
        }
        
        fxl <- ref[, names(para)[l]]
        gt <- matrix(0, nrow = n, ncol = nlam)
        gt[, id - offset] <- tmp * fxl
        g.alp2[[foo(j,l)]] <- t(gt %*% lam) %*% pr
      }
    }
  }
  
  if(length(g.alp2) == 0){
    g.alp2 <- NULL
  }
  
  g.alp2
  
}


