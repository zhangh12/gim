
# alp for each model: intercept, and other coefficients for covariates x (if any)
# by adding a constant column in data 'int', we do not distinguish intercept and other coefficients
# alp is NOT the 'alp' in code below, instead it is c(alp, bet)
# for each model, define 
# g = 
# (delta - delta_i) * phi(x)
# where delta = exp(X * theta)/(1 + exp(X * theta))
# delta_i = exp(X_1i * alp + X_2i * bet)/(exp(X_1i * alp + X_2i * bet)) for ith auxiliary model
gfunction.lo <- function(para, para.id, int){
  
  int$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.alp)
  
  the <- para[id.the$start[1]:id.the$end[1]]
  fx <- as.matrix(int[, names(the), drop = FALSE])
  
  n <- nrow(int)
  nt <- max(id.lam)
  g <- matrix(0, nrow = n, ncol = nt)
  offset <- min(id.alp)-1
  
  for(i in 1:nmodel){
    alp <- para[id.alp$start[i]:id.alp$end[i]]
    bet <- para[id.bet$start[i]:id.bet$end[i]]
    gam <- c(alp, bet)
    
    rx <- as.matrix(int[, names(gam), drop = FALSE])
    
    tmp1 <- as.vector(exp(fx %*% the))
    tmp2 <- as.vector(exp(rx %*% gam))
    delta <- tmp1/(1+tmp1) - tmp2/(1+tmp2)
    
    id <- c(id.alp$start[i]:id.alp$end[i], id.bet$start[i]:id.bet$end[i])
    g[, id - offset] <- rx[, names(para)[id], drop = FALSE] * delta
    
  }
  
  g
  
}
