

cov.cml <- function(para, map, data, ref, sample.info, V, bet0, outcome){
  
  para[map$all.bet] <- bet0
  names(para)[map$all.bet] <- names(bet0)
  Jv0 <- -hess.cc(para, map, data, ref, solve(V), bet0, sample.info, outcome)
  
  n <- nrow(data)
  n1 <- sum(data[, outcome])
  n0 <- n - n1
  
  Iv0 <- Jv0
  Iv0[] <- 0
  Iv0[map$lam, map$lam] <- -Jv0[map$lam, map$lam]
  Iv0[map$the, map$the] <- Jv0[map$the, map$the]
  Iv0[map$all.bet, map$all.bet] <- Jv0[map$all.bet, map$all.bet]
  v <- Jv0[, map$lam[1]]
  Iv0 <- Iv0 - 1/n * n1/n * n0/n * v %*% t(v)
  
  id <- c(map$lam, map$the, map$all.alp)
  mat <- Iv0[id, id] + Jv0[id, map$all.bet] %*% V %*% Jv0[map$all.bet, id]
  
  vcov0 <- solve(Jv0[id, id]) %*% mat %*% solve(Jv0[id, id])
  
  vcov0
  
}

