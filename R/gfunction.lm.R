
# alp for each model: tau (alp[1]), coefficients for covariates x (if any, including intercept)
# by adding a constant column in 'data', we do not distinguish intercept and other coefficients
# for each model, define 
# g = 
# (x * the - phi(x) * alp[-1])^2 + sigma - tau
# (x * the - phi(x) * alp[-1]) * phi(x)
gfunction.lm <- function(para, para.id, data){
  
  data$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.bet)
  
  sigma <- para[id.the$start[1]]
  the <- para[(id.the$start[1]+1):id.the$end[1]]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  
  n <- nrow(data)
  nlam <- max(id.lam)
  g <- matrix(0, nrow = n, ncol = nlam)
  offset <- max(id.the)
  
  for(i in 1:nmodel){
    id.tau <- id.alp$start[i]
    tau <- para[id.tau]
    
    id.a <- alp.index.lm(id.alp, i)
    alp.exist <- !is.null(id.a)
    if(alp.exist){
      alp <- para[id.a]
    }else{
      alp <- NULL
    }
    
    id.b <- id.bet$start[i]:id.bet$end[i]
    bet <- para[id.b]
    gam <- c(alp, bet)
    
    rx <- as.matrix(data[, names(gam), drop = FALSE])
    
    delta <- as.vector(fx %*% the - rx %*% gam)
    
    g[, id.tau - offset] <- delta^2 + sigma - tau
    id <- c(id.a, id.b)
    g[, id - offset] <- rx[, names(para)[id], drop = FALSE] * delta
    
    rm(id.tau, id.a, id.b, alp.exist)
  }
  
  g
  
}
