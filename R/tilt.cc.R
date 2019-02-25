
# compute Delta = exp(X * theta) for internal data
# compute delta_i = exp(X_1i * alp + X_2i * bet) for ith auxiliary model
tilt.cc <- function(para, map, ref){
  
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
  
  list(Delta = Delta, delta = delta)
  
}
