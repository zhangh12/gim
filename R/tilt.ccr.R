
# compute Delta = exp(X * theta) for external data
# compute delta_i = exp(X_1i * alp + X_2i * bet) for ith auxiliary model
tilt.ccr <- function(para, map, data, ref){
  
  the <- para[map$the]
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  Delta <- as.vector(exp(fx %*% the))
  
  nmodel <- length(map$bet)
  n <- nrow(ref)
  delta <- matrix(NA, nrow = n, ncol = nmodel)
  
  for(i in 1:nmodel){
    id <- c(alp.index.cc(map, i), map$bet[[i]])
    gam <- para[id]
    
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    delta[, i] <- as.vector(exp(rx %*% gam))
  }
  
  ome <- para[map$ome]
  fx0 <- as.matrix(data[, names(the), drop = FALSE])
  Delta0 <- as.vector(exp(ome + fx0 %*% the))
  
  list(Delta = Delta, delta = delta, Delta0 = Delta0)
  
}
