

gfunction.bet.xi.cc <- function(g.bet, xi){
  
  n <- nrow(g.bet[[1]])
  nbet <- length(g.bet)
  
  g.bet.xi <- matrix(NA, nrow = n, ncol = nbet)
  
  k <- 0
  for(m in g.bet){
    k <- k + 1
    g.bet.xi[, k] <- m[, -1, drop = FALSE] %*% xi
  }
  
  g.bet.xi
  
}
