

score.lm <- function(para, map, data, ref, inv.V, bet0, outcome){
  
  #sc1 <- grad(obj.lm, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome)
  #return(sc1)
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  n <- nrow(data)
  
  lam <- para[map$lam]
  
  sigma <- para[map$the[1]]
  the <- para[map$the[-1]]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  y <- data[, outcome]
  res <- y - as.vector(fx %*% the)
  
  g <- gfunction.lm(para, map, ref)
  
  pr <- as.vector(1/(1+g %*% lam))
  
  np <- length(para)
  
  sc <- rep(NA, np)
  names(sc) <- names(para)
  sc[map$lam] <- -as.vector(t(g) %*% pr)
  
  dlogL <- c(-n/2/sigma + 1/2/sigma^2 * sum(res^2), as.vector(t(fx) %*% res/sigma))
  g.the <- gfunction.the.lm(para, map, ref)
  id <- map$the
  for(i in 1:length(id)){
    tmp <- as.vector(g.the[[i]] %*% lam)
    sc[id[i]] <- dlogL[i] - sum(tmp * pr)
    rm(tmp)
  }
  
  g.alp <- gfunction.alp.lm(para, map, ref)
  k <- max(map$the)
  for(i in seq_along(g.alp)){
    tmp <- as.vector(g.alp[[i]] %*% lam)
    k <- k + 1
    sc[k] <- -sum(tmp * pr)
    rm(tmp)
  }
  
  bet <- para[map$all.bet]
  dqf <- as.vector(inv.V %*% (bet - bet0))
  g.bet <- gfunction.bet.lm(para, map, ref)
  k <- min(map$all.bet) - 1
  for(i in 1:length(g.bet)){
    tmp <- as.vector(g.bet[[i]] %*% lam)
    k <- k + 1
    sc[k] <- -sum(tmp * pr) - dqf[i]
  }
  
  sc
  
}

