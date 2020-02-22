

cov.ccr <- function(para, map, data, ref, sample.info, V, bet0, outcome){
  
  Jv0 <- -hess.ccr(para, map, data, ref, solve(V), bet0, sample.info, outcome)
  
  Iv0 <- Jv0
  Iv0[] <- 0
  Iv0[map$lam, map$lam] <- -Jv0[map$lam, map$lam]
  Iv0[-map$lam, -map$lam] <- Jv0[-map$lam, -map$lam]
  
  N <- nrow(ref)
  n <- nrow(data)
  n1 <- sum(data[, outcome])
  n0 <- n - n1
  rho <- n1 / n0
  v <- Jv0[, map$ome]
  Iv0 <- Iv0 - 1/n * (1 + rho)^2 / rho * v %*% t(v)
  vcov0 <- solve(Jv0) %*% Iv0 %*% solve(Jv0)
  
  vcov0
  
}

