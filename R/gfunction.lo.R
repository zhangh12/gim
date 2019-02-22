
# alp for each model: coefficients for covariates x (if any, including intercept)
# by adding a constant column in 'data', we do not distinguish intercept and other coefficients
# alp is NOT the 'alp' in code below, instead it is c(alp, bet)
# for each model, define 
# g = 
# (delta - delta_i) * phi(x)
# where delta = exp(X * theta)/(1 + exp(X * theta))
# delta_i = exp(X_1i * alp + X_2i * bet)/(exp(X_1i * alp + X_2i * bet)) for ith auxiliary model
gfunction.lo <- function(para, map, data){
  
  data$'(Intercept)' <- 1
  
  nmodel <- length(map$bet)
  
  the <- para[map$the]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  
  n <- nrow(data)
  nlam <- max(map$lam)
  g <- matrix(0, nrow = n, ncol = nlam)
  offset <- max(map$the)
  for(i in 1:nmodel){
    
    id.a <- alp.index.lo(map, i)
    alp.exist <- !is.null(id.a)
    if(alp.exist){
      alp <- para[id.a]
    }else{
      alp <- NULL
    }
    
    id.b <- map$bet[[i]]
    bet <- para[id.b]
    gam <- c(alp, bet)
    
    rx <- as.matrix(data[, names(gam), drop = FALSE])
    
    tmp1 <- as.vector(exp(fx %*% the))
    tmp2 <- as.vector(exp(rx %*% gam))
    delta <- tmp1/(1+tmp1) - tmp2/(1+tmp2)
    
    id <- c(id.a, id.b)
    g[, id - offset] <- rx[, names(para)[id], drop = FALSE] * delta
    
    rm(id.a, id.b, alp.exist)
    
  }
  
  g
  
}

