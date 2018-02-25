

gfunction.bet.lm <- function(para, para.id, int){
  
  int$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.alp)
  
  sigma <- para[id.the$start[1]]
  the <- para[(id.the$start[1]+1):id.the$end[1]]
  fx <- as.matrix(int[, names(the), drop = FALSE])
  
  g.bet <- list()
  
  n <- nrow(int)
  nt <- max(id.lam)
  offset <- min(id.alp)-1
  
  k <- 0
  for(i in 1:nmodel){
    
    tau <- para[id.alp$start[i]]
    alp <- para[(id.alp$start[i] + 1):id.alp$end[i]]
    bet <- para[id.bet$start[i]:id.bet$end[i]]
    gam <- c(alp, bet)
    
    rx <- as.matrix(int[, names(gam), drop = FALSE])
    
    delta <- as.vector(fx %*% the - rx %*% gam)
    
    id <- id.bet$start[i]:id.bet$end[i]
    for(j in id){
      rx0 <- rx[, names(para)[j]]
      gb <- matrix(0, nrow = n, ncol = nt)
      gb[, id.alp$start[i] - offset] <- -2 * delta * rx0
      gb[, (id.alp$start[i]+1):(id.alp$end[i]) - offset] <- -rx[, names(alp), drop = FALSE] * rx0
      gb[, id - offset] <- -rx[, names(bet), drop = FALSE] * rx0
      k <- k + 1
      g.bet[[k]] <- gb
      rm(gb)
    }
  }
  
  g.bet
  
}

