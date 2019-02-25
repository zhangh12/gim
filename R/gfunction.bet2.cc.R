

gfunction.bet2.cc <- function(para, map, ref, Delta, delta, ncase, nctrl){
  
  nmodel <- length(map$bet)
  
  g.bet2 <- list()
  
  ref <- as.matrix(ref)
  n <- nrow(ref)
  
  const <- list()
  for(i in 1:nmodel){
    id <- c(alp.index.cc(map, i), map$bet[[i]])
    gam <- para[id]
    rx <- ref[, names(gam), drop = FALSE]
    
    rho.i <- ncase[i, i] / nctrl[i, i]
    const[[i]] <- -rx * (delta[, i] * (1 - rho.i * delta[, i]) * (1 + rho.i * Delta) / (1 + rho.i * delta[, i])^3)
  }
  
  nlam <- max(map$lam)
  offset <- max(map$the) - 1
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  for(i in 1:nmodel){
    id.a <- alp.index.cc(map, i)
    id.b <- map$bet[[i]]
    
    id <- c(id.a, id.b)
    for(j in id.b){
      fxj <- ref[, names(para)[j]]
      tmp <- const[[i]] * fxj
      for(l in id.b){
        if(l < j){
          next
        }
        
        fxl <- ref[, names(para)[l]]
        gt <- matrix(0, nrow = n, ncol = nlam)
        gt[, id - offset] <- tmp * fxl
        g.bet2[[foo(j,l)]] <- gt
      }
    }
  }
  
  g.bet2
  
}

