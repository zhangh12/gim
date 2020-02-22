
score.ccr <- function(para, map, data, ref, inv.V, bet0, sample.info, outcome){
  
  #return(grad(obj.ccr, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome))
  
  ncase <- sample.info$ncase
  nctrl <- sample.info$nctrl
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  
  lam <- para[map$lam]
  xi <- lam[-1]
  
  ome <- para[map$ome]
  the <- para[map$the]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  y <- data[, outcome]
  n1 <- sum(y)
  n0 <- sum(!y)
  rho <- n1/n0
  
  tilt <- tilt.ccr(para, map, data, ref)
  Delta <- tilt$Delta
  delta <- tilt$delta
  Delta0 <- tilt$Delta0
  
  g <- gfunction.cc(para, map, ref, Delta, delta, ncase, nctrl)
  
  pr <- as.vector(1/(1+g %*% lam))
  
  np <- length(para)
  
  sc <- rep(NA, np)
  names(sc) <- names(para)
  sc[map$lam] <- -as.vector(t(g) %*% pr)
  
  pred <- 1 - 1/(1 + rho * Delta0)
  sc[map$ome] <- sum(y - pred)
  
  dlogL <- as.vector(t(fx) %*% (y - pred))
  g.the <- gfunction.the.cc(para, map, ref, Delta, delta, ncase, nctrl)
  g.the.xi <- gfunction.the.xi.cc(g.the, xi)
  tmp <- g.the.xi + lam[1] * Delta * ref[, names(the)]
  sc[map$the] <- dlogL - t(tmp) %*% pr
  rm(tmp)
  
  if(!is.null(map$all.alp)){
    g.alp <- gfunction.alp.cc(para, map, ref, Delta, delta, ncase, nctrl)
    g.alp.xi <- gfunction.alp.xi.cc(g.alp, xi)
    sc[map$all.alp] <- -t(g.alp.xi) %*% pr
  }
  
  bet <- para[map$all.bet]
  dqf <- as.vector(inv.V %*% (bet - bet0))
  g.bet <- gfunction.bet.cc(para, map, ref, Delta, delta, ncase, nctrl)
  g.bet.xi <- gfunction.alp.xi.cc(g.bet, xi)
  sc[map$all.bet] <- -t(g.bet.xi) %*% pr - dqf
  
  sc
  
}

