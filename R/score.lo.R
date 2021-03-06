

score.lo <- function(para, map, data, ref, inv.V, bet0, outcome){
  
  #return(grad(obj.lo, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome))
  
  data$'(Intercept)' <- 1
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  n <- nrow(data)
  
  lam <- para[map$lam]
  
  the <- para[map$the]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  y <- data[, outcome]
  lin <- as.vector(fx %*% the)
  elin <- exp(lin)
  res <- y-elin/(1+elin)
  
  g <- gfunction.lo(para, map, ref)
  
  pr <- as.vector(1/(1+g %*% lam))
  
  np <- length(para)
  
  sc <- rep(NA, np)
  names(sc) <- names(para)
  sc[map$lam] <- -as.vector(t(g) %*% pr)
  
  dlogL <- as.vector(t(fx) %*% res)
  g.the <- gfunction.the.lo(para, map, ref)
  id <- map$the
  for(i in 1:length(id)){
    tmp <- as.vector(g.the[[i]] %*% lam)
    sc[id[i]] <- dlogL[i] - sum(tmp * pr)
    rm(tmp)
  }
  
  g.alp <- gfunction.alp.lo(para, map, ref)
  k <- max(map$the)
  for(ga in g.alp){
    tmp <- as.vector(ga %*% lam)
    k <- k + 1
    sc[k] <- -sum(tmp * pr)
    rm(tmp)
  }
  
  bet <- para[map$all.bet]
  dqf <- as.vector(inv.V %*% (bet - bet0))
  g.bet <- gfunction.bet.lo(para, map, ref)
  k <- min(map$all.bet) - 1
  for(i in 1:length(g.bet)){
    tmp <- as.vector(g.bet[[i]] %*% lam)
    k <- k + 1
    sc[k] <- -sum(tmp * pr) - dqf[i]
  }
  
  sc
  
}

