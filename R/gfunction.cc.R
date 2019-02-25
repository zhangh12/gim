
# alp for each model: coefficients for covariates x (if any, including intercept)
# by adding a constant column in 'data', we do not distinguish intercept and other coefficients
# alp is NOT the 'alp' in code below, instead it is c(alp, bet)
# for each model, define 
# g = 
# (Delta - delta_i) / (1 + rho_i * delta_i) X_i
# where Delta = exp(X * theta)
# delta_i = exp(X_1i * alp + X_2i * bet) for ith auxiliary model
# final g = [Delta - 1, g]
gfunction.cc <- function(para, map, ref, Delta, delta, ncase, nctrl){
  
  nmodel <- length(map$bet)
  
  the <- para[map$the]
  
  n <- nrow(ref)
  nlam <- max(map$lam)
  g <- matrix(0, nrow = n, ncol = nlam)
  
  g[, 1] <- Delta - 1
  offset <- max(map$the) - 1
  for(i in 1:nmodel){
    
    id <- c(alp.index.cc(map, i), map$bet[[i]])
    gam <- para[id]
    
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    
    rho.i <- ncase[i, i] / nctrl[i, i]
    
    g[, id - offset] <- rx * ((Delta - delta[, i]) / (1 + rho.i * delta[, i]))
    
  }
  
  g
  
}

