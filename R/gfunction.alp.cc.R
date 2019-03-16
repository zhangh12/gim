

gfunction.alp.cc <- function(para, map, ref, Delta, delta, ncase, nctrl){
  
  nmodel <- length(map$bet)
  
  g.alp <- list()
  
  n <- nrow(ref)
  nlam <- max(map$lam)
  offset <- max(map$the)
  
  k <- 0
  for(i in 1:nmodel){
    
    id.a <- alp.index.cc(map, i)
    alp.exist <- !is.null(id.a)
    if(alp.exist){
      alp <- para[id.a]
    }else{
      #alp <- NULL
      next
    }
    
    id.b <- map$bet[[i]]
    id <- c(id.a, id.b)
    gam <- para[id]
    
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    rho.i <- ncase[i, i] / nctrl[i, i]
    const <- rx * (delta[, i] * (1 + rho.i * Delta) / (1 + rho.i * delta[, i])^2)
    
    for(j in id.a){
      rx0 <- rx[, names(para)[j]]
      ga <- matrix(0, nrow = n, ncol = nlam - 1)
      ga[, id - offset] <- -const * rx0
      k <- k + 1
      g.alp[[k]] <- ga
      rm(ga)
    }
  }
  
  if(length(g.alp) == 0){
    g.alp <- NULL
  }
  
  g.alp
  
}

