

cov.lo <- function(para, map, data, ref, model, nsample, V, bet0, outcome){
  
  data$'(Intercept)' <- 1
  
  g <- gfunction.lo(para, map, ref)
  g.the <- gfunction.the.lo(para, map, ref)
  g.alp <- gfunction.alp.lo(para, map, ref)
  g.bet <- gfunction.bet.lo(para, map, ref)
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  n <- nrow(data)
  m <- nrow(ref)
  r <- m/n
  
  lam <- para[map$lam]
  pr <- as.vector(1/(1+g %*% lam))
  pr <- pr/sum(pr)
  pr0 <- 1/n
  
  J.tt <- -(t(g) %*% (g * pr)) * r
  
  nthe <- length(g.the)
  J.tthe <- matrix(0, nrow = nlam, ncol = nthe)
  for(i in 1:nthe){
    J.tthe[, i] <- t(g.the[[i]]) %*% pr * r
  }
  
  nalp <- length(g.alp)
  if(nalp > 0){
    J.talp <- matrix(0, nrow = nlam, ncol = nalp)
    for(i in 1:nalp){
      J.talp[, i] <- t(g.alp[[i]]) %*% pr * r
    }
  }
  
  nbet <- length(g.bet)
  J.tbet <- matrix(0, nrow = nlam, ncol = nbet)
  for(i in 1:nbet){
    J.tbet[, i] <- t(g.bet[[i]]) %*% pr * r
  }
  
  the <- para[map$the]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  elin <- as.vector(exp(fx %*% the))
  yh <- elin/(1+elin)
  J.the <- t(fx) %*% (fx * (yh * (1-yh) * pr0))
  
  suppressMessages(Sigma0 <- Sigma0.lo(para, map, ref, model, nsample, outcome))
  J.bet <- solve(Sigma0)/n
  #J.bet <- solve(V)/n
  
  np <- length(para)
  Jv <- matrix(0, nrow = np, ncol = np)
  Iv <- matrix(0, nrow = np, ncol = np)
  
  Jv[map$lam, map$lam] <- J.tt
  Jv[map$lam, map$the] <- J.tthe
  Jv[map$the, map$lam] <- t(J.tthe)
  if(nalp > 0){
    Jv[map$lam, map$all.alp] <- J.talp
    Jv[map$all.alp, map$lam] <- t(J.talp)
  }
  Jv[map$lam, map$all.bet] <- J.tbet
  Jv[map$all.bet, map$lam] <- t(J.tbet)
  Jv[map$the, map$the] <- J.the
  Jv[map$all.bet, map$all.bet] <- J.bet
  
  Iv[map$lam, map$lam] <- -J.tt
  Iv[map$the, map$the] <- J.the
  Iv[map$all.bet, map$all.bet] <- J.bet
  #Iv[map$all.bet, map$all.bet] <- n * J.bet %*% Sigma0 %*% J.bet
  
  vcov <- solve(Jv) %*% Iv %*% solve(Jv)/n
  
  vcov
  
  #################
  # seems like this sample version is better
  # so abandent everything before this line
  #################
  
  Jv0 <- -hess.lo(para, map, data, ref, solve(V), bet0, outcome)/n
  
  Iv0 <- matrix(0, nrow = np, ncol = np)
  Iv0[map$lam, map$lam] <- -Jv0[map$lam, map$lam]
  Iv0[map$the, map$the] <- Jv0[map$the, map$the]
  Iv0[map$all.bet, map$all.bet] <- Jv0[map$all.bet, map$all.bet]
  vcov0 <- solve(Jv0) %*% Iv0 %*% solve(Jv0)/n
  
  vcov0
  
}

