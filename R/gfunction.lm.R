
# alp for each model: tau (alp[1]), intercept, and other coefficients for covariates x (if any)
# by adding a constant column in data 'int', we do not distinguish intercept and other coefficients
# for each model, define 
# g = 
# (x * the - phi(x) * alp[-1])^2 + sigma - tau
# (x * the - phi(x) * alp[-1]) * phi(x)
gfunction.lm <- function(para, para.id, int){
  
  int$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.alp)
  
  sigma <- para[id.the$start[1]]
  the <- para[(id.the$start[1]+1):id.the$end[1]]
  fx <- as.matrix(int[, names(the), drop = FALSE])
  
  n <- nrow(int)
  nt <- max(id.lam)
  g <- matrix(0, nrow = n, ncol = nt)
  offset <- min(id.alp)-1
  
  for(i in 1:nmodel){
    tau <- para[id.alp$start[i]]
    alp <- para[(id.alp$start[i] + 1):id.alp$end[i]]
    bet <- para[id.bet$start[i]:id.bet$end[i]]
    gam <- c(alp, bet)
    
    rx <- as.matrix(int[, names(gam), drop = FALSE])
    
    delta <- as.vector(fx %*% the - rx %*% gam)
    
    g[, id.alp$start[i] - offset] <- delta^2 + sigma - tau
    id <- c((id.alp$start[i]+1):id.alp$end[i], id.bet$start[i]:id.bet$end[i])
    g[, id - offset] <- rx[, names(para)[id], drop = FALSE] * delta
    
  }
  
  g
  
}
