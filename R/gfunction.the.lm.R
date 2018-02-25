

gfunction.the.lm <- function(para, para.id, int){
  
  int$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.alp)
  
  sigma <- para[id.the$start[1]]
  the <- para[(id.the$start[1]+1):id.the$end[1]]
  fx <- as.matrix(int[, names(the), drop = FALSE])
  
  g.the <- list()
  
  nthe <- length(the)
  n <- nrow(int)
  
  nt <- max(id.lam)
  
  gt <- matrix(0, nrow = n, ncol = nt)
  offset <- min(id.alp)-1
  for(i in 1:nmodel){
    gt[, id.alp$start[i] - offset] <- 1
  }
  
  g.the[[1]] <- gt
  rm(gt)
  
  for(j in 1:nthe){
    gt <- matrix(0, nrow = n, ncol = nt)
    fx0 <- fx[, names(the)[j]]
    for(i in 1:nmodel){
      tau <- para[id.alp$start[i]]
      alp <- para[(id.alp$start[i] + 1):id.alp$end[i]]
      bet <- para[id.bet$start[i]:id.bet$end[i]]
      gam <- c(alp, bet)
      
      rx <- as.matrix(int[, names(gam), drop = FALSE])
      
      delta <- as.vector(fx %*% the - rx %*% gam)
      
      gt[, id.alp$start[i] - offset] <- 2 * delta * fx0
      id <- c((id.alp$start[i]+1):id.alp$end[i], id.bet$start[i]:id.bet$end[i])
      gt[, id - offset] <- rx[, names(para)[id], drop = FALSE] * fx0
      
    }
    g.the[[j+1]] <- gt
    rm(gt)
  }
  
  g.the
  
}

