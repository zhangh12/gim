

gfunction.bet.cc <- function(para, map, ref, Delta, delta, ncase, nctrl){
  
  nmodel <- length(map$bet)
  
  g.bet <- list()
  
  n <- nrow(ref)
  nlam <- max(map$lam)
  offset <- max(map$the) - 1
  
  k <- 0
  for(i in 1:nmodel){
    
    id.a <- alp.index.lo(map, i)
    alp.exist <- !is.null(id.a)
    if(alp.exist){
      alp <- para[id.a]
    }else{
      alp <- NULL
    }
    
    id.b <- map$bet[[i]]
    id <- c(id.a, id.b)
    gam <- para[id]
    
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    rho.i <- ncase[i, i] / nctrl[i, i]
    const <- rx * (delta[, i] * (1 + rho.i * Delta) / (1 + rho.i * delta[, i])^2)
    
    for(j in id.b){
      rx0 <- rx[, names(para)[j]]
      gb <- matrix(0, nrow = n, ncol = nlam)
      gb[, id - offset] <- -const * rx0
      k <- k + 1
      g.bet[[k]] <- gb
      rm(gb)
    }
  }
  
  g.bet
  
}


