

gfunction.bet.lam.lo <- function(g.bet, lam){
  
  n <- nrow(g.bet[[1]])
  nbet <- length(g.bet)
  g.bet.lam <- matrix(NA, nrow = n, ncol = nbet)
  k <- 0
  for(m in g.bet){
    k <- k + 1
    g.bet.lam[, k] <- m %*% lam
  }
  
  g.bet.lam
  
}
