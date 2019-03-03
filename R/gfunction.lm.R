
# alp for each model: tau (alp[1]), coefficients for covariates x (if any, including intercept)
# by adding a constant column in 'ref', we do not distinguish intercept and other coefficients
# for each model, define 
# g = 
# (x * the - phi(x) * alp[-1])^2 + sigma - tau
# (x * the - phi(x) * alp[-1]) * phi(x)
gfunction.lm <- function(para, map, ref){
  
  ref$'(Intercept)' <- 1
  
  nmodel <- length(map$bet)
  
  sigma <- para[map$the[1]]
  the <- para[map$the[-1]]
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  
  n <- nrow(ref)
  nlam <- max(map$lam)
  g <- matrix(0, nrow = n, ncol = nlam)
  offset <- max(map$the)
  for(i in 1:nmodel){
    id.tau <- map$alp[[i]][1]
    tau <- para[id.tau]
    
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
    
    g[, id.tau - offset] <- delta^2 + sigma - tau
    id <- c(id.a, id.b)
    g[, id - offset] <- rx[, names(para)[id], drop = FALSE] * delta
    
    rm(id.tau, id.a, id.b, alp.exist)
  }
  
  g
  
}

