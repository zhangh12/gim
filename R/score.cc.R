
score.cc <- function(para, map, data, ref, inv.V, bet0, sample.info, outcome){
  
  #return(grad(obj.cc, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome))
  
  ncase <- sample.info$ncase
  nctrl <- sample.info$nctrl
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  
  lam <- para[map$lam]
  
  the <- para[map$the]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  y <- data[, outcome]
  
  tilt <- tilt.cc(para, map, ref)
  Delta <- tilt$Delta
  delta <- tilt$delta
  
  g <- gfunction.cc(para, map, ref, Delta, delta, ncase, nctrl)
  
  pr <- as.vector(1/(1+g %*% lam))
  
  np <- length(para)
  
  sc <- rep(NA, np)
  names(sc) <- names(para)
  sc[map$lam] <- -as.vector(t(g) %*% pr)
  
  dlogL <- as.vector(t(fx) %*% y)
  g.the <- gfunction.the.cc(para, map, ref, Delta, delta, ncase, nctrl)
  id <- map$the
  for(i in 1:length(id)){
    tmp <- as.vector(g.the[[i]] %*% lam)
    sc[id[i]] <- dlogL[i] - sum(tmp * pr)
    rm(tmp)
  }
  
  g.alp <- gfunction.alp.cc(para, map, ref, Delta, delta, ncase, nctrl)
  k <- max(map$the)
  for(ga in g.alp){
    tmp <- as.vector(ga %*% lam)
    k <- k + 1
    sc[k] <- -sum(tmp * pr)
    rm(tmp)
  }
  
  bet <- para[map$all.bet]
  dqf <- as.vector(inv.V %*% (bet - bet0))
  g.bet <- gfunction.bet.cc(para, map, ref, Delta, delta, ncase, nctrl)
  k <- min(map$all.bet) - 1
  for(i in 1:length(g.bet)){
    tmp <- as.vector(g.bet[[i]] %*% lam)
    k <- k + 1
    sc[k] <- -sum(tmp * pr) - dqf[i]
  }
  
  sc
  
}

