

gfunction.alp.lo <- function(para, para.id, int){
  
  int$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.alp)
  
  the <- para[id.the$start[1]:id.the$end[1]]
  fx <- as.matrix(int[, names(the), drop = FALSE])
  
  g.alp <- list()
  
  n <- nrow(int)
  nt <- max(id.lam)
  offset <- min(id.alp)-1
  
  k <- 0
  for(i in 1:nmodel){
    
    alp <- para[id.alp$start[i]:id.alp$end[i]]
    bet <- para[id.bet$start[i]:id.bet$end[i]]
    gam <- c(alp, bet)
    
    rx <- as.matrix(int[, names(gam), drop = FALSE])
    tmp <- as.vector(exp(rx %*% gam))
    yi <- tmp/(1+tmp)
    
    id <- id.alp$start[i]:id.alp$end[i]
    for(j in id){
      rx0 <- rx[, names(para)[j]]
      ga <- matrix(0, nrow = n, ncol = nt)
      
      ga[, id - offset] <- -rx[, names(alp), drop = FALSE] * (rx0 * yi * (1-yi))
      ga[, (id.bet$start[i]):id.bet$end[i] - offset] <- -rx[, names(bet), drop = FALSE] * (rx0 * yi * (1-yi))
      k <- k + 1
      g.alp[[k]] <- ga
      rm(ga)
    }
  }
  
  g.alp
  
}

