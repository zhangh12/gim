

cov.lm <- function(para, map, data, ref, model, nsample, V, bet0, outcome){
  
  #################
  # seems like this sample version is better
  #################
  
  n <- nrow(data)
  np <- length(para)
  Jv0 <- -hess.lm(para, map, data, ref, solve(V), bet0, outcome)/n
  
  Iv0 <- matrix(0, nrow = np, ncol = np)
  Iv0[map$lam, map$nlam] <- -Jv0[map$lam, map$lam]
  Iv0[map$the, map$the] <- Jv0[map$the, map$the]
  Iv0[map$all.bet, map$all.bet] <- Jv0[map$all.bet, map$all.bet]
  vcov0 <- solve(Jv0) %*% Iv0 %*% solve(Jv0)/n
  
  vcov0
  
}

