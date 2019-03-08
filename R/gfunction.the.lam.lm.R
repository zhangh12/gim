
# the first component in g.the.lam should be all zero
gfunction.the.lam.lm <- function(g.the, lam){
  
  n <- nrow(g.the[[1]])
  nthe <- length(g.the)
  g.the.lam <- matrix(NA, nrow = n, ncol = nthe)
  k <- 0
  for(m in g.the){
    k <- k + 1
    g.the.lam[, k] <- m %*% lam
  }
  
  g.the.lam
  
}
