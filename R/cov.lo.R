

cov.lo <- function(para, para.id, data, model, nsample, outcome = 'y'){
  
  data$'(Intercept)' <- 1
  
  g <- gfunction.lo(para, para.id, data)
  g.the <- gfunction.the.lo(para, para.id, data)
  g.alp <- gfunction.alp.lo(para, para.id, data)
  g.bet <- gfunction.bet.lo(para, para.id, data)
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.bet)
  nlam <- max(id.lam)
  n <- nrow(data)
  
  lam <- para[id.lam$start[1]:id.lam$end[1]]
  pr <- as.vector(1/(1+g %*% lam))
  pr <- pr/sum(pr)
  
  J.tt <- -(t(g) %*% (g * pr))
  
  nthe <- length(g.the)
  J.tthe <- matrix(0, nrow = nlam, ncol = nthe)
  for(i in 1:nthe){
    J.tthe[, i] <- t(g.the[[i]]) %*% pr
  }
  
  nalp <- length(g.alp)
  if(nalp > 0){
    J.talp <- matrix(0, nrow = nlam, ncol = nalp)
    for(i in 1:nalp){
      J.talp[, i] <- t(g.alp[[i]]) %*% pr
    }
  }
  
  nbet <- length(g.bet)
  J.tbet <- matrix(0, nrow = nlam, ncol = nbet)
  for(i in 1:nbet){
    J.tbet[, i] <- t(g.bet[[i]]) %*% pr
  }
  
  the <- para[id.the$start[1]:id.the$end[1]]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  elin <- as.vector(exp(fx %*% the))
  yh <- elin/(1+elin)
  J.the <- t(fx) %*% (fx * (yh * (1-yh) * pr))
  
  suppressMessages(Sigma0 <- Sigma0.lo(para, para.id, data, model, nsample, outcome))
  J.bet <- solve(Sigma0)/n
  
  np <- length(para)
  Jv <- matrix(0, nrow = np, ncol = np)
  Iv <- matrix(0, nrow = np, ncol = np)
  
  Jv[1:nlam, 1:nlam] <- J.tt
  Jv[1:nlam, min(id.the):max(id.the)] <- J.tthe
  Jv[min(id.the):max(id.the), 1:nlam] <- t(J.tthe)
  if(nalp > 0){
    Jv[1:nlam, min(id.alp, na.rm = TRUE):max(id.alp, na.rm = TRUE)] <- J.talp
    Jv[min(id.alp, na.rm = TRUE):max(id.alp, na.rm = TRUE), 1:nlam] <- t(J.talp)
  }
  Jv[1:nlam, min(id.bet):max(id.bet)] <- J.tbet
  Jv[min(id.bet):max(id.bet), 1:nlam] <- t(J.tbet)
  Jv[min(id.the):max(id.the), min(id.the):max(id.the)] <- J.the
  Jv[min(id.bet):max(id.bet), min(id.bet):max(id.bet)] <- J.bet
  
  Iv[1:nlam, 1:nlam] <- -J.tt
  Iv[min(id.the):max(id.the), min(id.the):max(id.the)] <- J.the
  Iv[min(id.bet):max(id.bet), min(id.bet):max(id.bet)] <- J.bet
  
  vcov <- solve(Jv) %*% Iv %*% solve(Jv)/n
  
  vcov
  
}

