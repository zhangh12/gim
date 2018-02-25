

cov.lm <- function(para, para.id, int, model, nsample){
  
  int$'(Intercept)' <- 1
  
  g <- gfunction.lm(para, para.id, int)
  g.the <- gfunction.the.lm(para, para.id, int)
  g.alp <- gfunction.alp.lm(para, para.id, int)
  g.bet <- gfunction.bet.lm(para, para.id, int)
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.alp)
  nt <- max(id.lam)
  n <- nrow(int)
  
  lam <- para[id.lam$start[1]:id.lam$end[1]]
  pr <- as.vector(1/(1+g %*% lam))
  pr <- pr/sum(pr)
  
  J.tt <- -(t(g) %*% (g * pr))
  
  nthe <- length(g.the)
  J.tthe <- matrix(0, nrow = nt, ncol = nthe)
  for(i in 1:nthe){
    J.tthe[, i] <- t(g.the[[i]]) %*% pr
  }
  
  nalp <- length(g.alp)
  J.talp <- matrix(0, nrow = nt, ncol = nalp)
  for(i in 1:nalp){
    J.talp[, i] <- t(g.alp[[i]]) %*% pr
  }
  
  nbet <- length(g.bet)
  J.tbet <- matrix(0, nrow = nt, ncol = nbet)
  for(i in 1:nbet){
    J.tbet[, i] <- t(g.bet[[i]]) %*% pr
  }
  
  J.the <- matrix(0, nrow = nthe, ncol = nthe)
  sigma <- para[id.the$start[1]]
  the <- para[(id.the$start[1]+1):id.the$end[1]]
  fx <- as.matrix(int[, names(the), drop = FALSE])
  J.the[1, 1] <- 1/2/sigma
  J.the[-1, -1] <- t(fx) %*% (fx * pr)
  J.the <- 1/sigma * J.the
  
  suppressMessages(Sigma0 <- Sigma0.lm(para, para.id, int, model, nsample))
  J.bet <- solve(Sigma0)/n
  
  np <- length(para)
  Jv <- matrix(0, nrow = np, ncol = np)
  Iv <- matrix(0, nrow = np, ncol = np)
  
  Jv[1:nt, 1:nt] <- J.tt
  Jv[1:nt, min(id.the):max(id.the)] <- J.tthe
  Jv[min(id.the):max(id.the), 1:nt] <- t(J.tthe)
  Jv[1:nt, min(id.alp):max(id.alp)] <- J.talp
  Jv[min(id.alp):max(id.alp), 1:nt] <- t(J.talp)
  Jv[1:nt, min(id.bet):max(id.bet)] <- J.tbet
  Jv[min(id.bet):max(id.bet), 1:nt] <- t(J.tbet)
  Jv[min(id.the):max(id.the), min(id.the):max(id.the)] <- J.the
  Jv[min(id.bet):max(id.bet), min(id.bet):max(id.bet)] <- J.bet
  
  Iv[1:nt, 1:nt] <- -J.tt
  Iv[min(id.the):max(id.the), min(id.the):max(id.the)] <- J.the
  Iv[min(id.bet):max(id.bet), min(id.bet):max(id.bet)] <- J.bet
  
  vcov <- solve(Jv) %*% Iv %*% solve(Jv)/n
  
  vcov
  
}

