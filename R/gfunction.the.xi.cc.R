
gfunction.the.xi.cc <- function(g.the, xi){
  
  n <- nrow(g.the[[1]])
  nthe <- length(g.the)
  g.the.xi <- matrix(NA, nrow = n, ncol = nthe)
  k <- 0
  for(m in g.the){
    k <- k + 1
    #g.the.xi[, k] <- m[, -1, drop = FALSE] %*% xi
    g.the.xi[, k] <- m %*% xi
  }
  
  g.the.xi
  
}
