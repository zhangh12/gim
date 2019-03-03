

gfunction.bet.lam.lm <- function(g.bet, lam){
  
  n <- nrow(g.bet[[1]])
  nalp <- length(g.bet)
  g.bet.lam <- matrix(NA, nrow = n, ncol = nalp)
  k <- 0
  for(m in g.bet){
    k <- k + 1
    g.bet.lam[, k] <- m %*% lam
  }
  
  g.bet.lam
  
}
