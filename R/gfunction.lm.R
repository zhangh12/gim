
# alp for each model: coefficients for covariates x (if any, including intercept), 
# from v0.22.0, no tau in alp anymore, we do not estimate it
# by adding a constant column in 'ref', we do not distinguish intercept and other coefficients
# for each model, define 
# g = 
# (x * the - phi(x) * alp) * phi(x)
gfunction.lm <- function(para, map, ref){
  
  nmodel <- length(map$bet)
  
  sigma <- para[map$the[1]]
  the <- para[map$the[-1]]
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  
  n <- nrow(ref)
  nlam <- max(map$lam)
  g <- matrix(0, nrow = n, ncol = nlam)
  offset <- max(map$the)
  for(i in 1:nmodel){
    id.a <- alp.index.lm(map, i)
    alp.exist <- !is.null(id.a)
    if(alp.exist){
      alp <- para[id.a]
    }else{
      alp <- NULL
    }
    
    id.b <- map$bet[[i]]
    bet <- para[id.b]
    gam <- c(alp, bet)
    
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    
    delta <- as.vector(fx %*% the - rx %*% gam)
    
    id <- c(id.a, id.b)
    g[, id - offset] <- rx[, names(para)[id], drop = FALSE] * delta
    
  }
  
  g
  
}

