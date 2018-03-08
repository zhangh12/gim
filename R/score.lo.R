

score.lo <- function(para, para.id, data, inv.V, bet0, outcome){
  
  #return(grad(obj.lo, para, para.id = para.id, data = data, inv.V = inv.V, bet0 = bet0, outcome = outcome))
  
  data$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.bet)
  nlam <- max(id.lam)
  n <- nrow(data)
  
  lam <- para[id.lam$start[1]:id.lam$end[1]]
  
  the <- para[id.the$start[1]:id.the$end[1]]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  y <- data[, outcome]
  lin <- as.vector(fx %*% the)
  elin <- exp(lin)
  res <- y-elin/(1+elin)
  
  g <- gfunction.lo(para, para.id, data)
  
  pr <- as.vector(1/(1+g %*% lam))
  
  np <- length(para)
  
  sc <- rep(NA, np)
  names(sc) <- names(para)
  sc[id.lam$start[1]:id.lam$end[1]] <- -as.vector(t(g) %*% pr)
  
  dlogL <- as.vector(t(fx) %*% res)
  g.the <- gfunction.the.lo(para, para.id, data)
  id <- id.the$start[1]:id.the$end[1]
  for(i in 1:length(id)){
    tmp <- as.vector(g.the[[i]] %*% lam)
    sc[id[i]] <- dlogL[i] - sum(tmp * pr)
    rm(tmp)
  }
  
  g.alp <- gfunction.alp.lo(para, para.id, data)
  k <- max(id.the)
  for(ga in g.alp){
    tmp <- as.vector(ga %*% lam)
    k <- k + 1
    sc[k] <- -sum(tmp * pr)
    rm(tmp)
  }
  
  bet <- para[min(id.bet):max(id.bet)]
  dqf <- as.vector(inv.V %*% (bet - bet0))
  g.bet <- gfunction.bet.lo(para, para.id, data)
  k <- min(id.bet) - 1
  for(i in 1:length(g.bet)){
    tmp <- as.vector(g.bet[[i]] %*% lam)
    k <- k + 1
    sc[k] <- -sum(tmp * pr) - dqf[i]
  }
  
  sc
  
}

