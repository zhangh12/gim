

gfunction.tau.lam.lm <- function(g.tau, lam){
  
  n <- nrow(g.tau[[1]])
  ntau <- length(g.tau)
  g.tau.lam <- matrix(NA, nrow = n, ncol = ntau)
  k <- 0
  for(m in g.tau){
    k <- k + 1
    id <- which(m[1, ] == -1)
    g.tau.lam[, k] <- m[, id] * lam[id]
  }
  
  g.tau.lam
  
}
