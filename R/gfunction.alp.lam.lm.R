

gfunction.alp.lam.lm <- function(g.alp, lam){
  
  n <- nrow(g.alp[[1]])
  nalp <- length(g.alp)
  g.alp.lam <- matrix(NA, nrow = n, ncol = nalp)
  k <- 0
  for(m in g.alp){
    k <- k + 1
    g.alp.lam[, k] <- m %*% lam
  }
  
  g.alp.lam
  
}
