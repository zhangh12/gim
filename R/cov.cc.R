

cov.cc <- function(para, map, data, ref, sample.info, V, bet0, outcome){
  
  Jv0 <- -hess.cc(para, map, data, ref, solve(V), bet0, sample.info, outcome)
  
  Iv0 <- Jv0
  Iv0[] <- 0
  Iv0[map$lam, map$lam] <- -Jv0[map$lam, map$lam]
  Iv0[map$the, map$the] <- Jv0[map$the, map$the]
  Iv0[map$all.bet, map$all.bet] <- Jv0[map$all.bet, map$all.bet]
  vcov0 <- solve(Jv0) %*% Iv0 %*% solve(Jv0)
  
  vcov0
  
}

