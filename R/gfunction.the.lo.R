

gfunction.the.lo <- function(para, para.id, int){
  
  int$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.alp)
  
  the <- para[id.the$start[1]:id.the$end[1]]
  fx <- as.matrix(int[, names(the), drop = FALSE])
  
  g.the <- list()
  
  nthe <- length(the)
  n <- nrow(int)
  
  nt <- max(id.lam)
  offset <- min(id.alp)-1
  
  for(j in 1:nthe){
    gt <- matrix(0, nrow = n, ncol = nt)
    fx0 <- fx[, names(the)[j]]
    for(i in 1:nmodel){
      alp <- para[id.alp$start[i]:id.alp$end[i]]
      bet <- para[id.bet$start[i]:id.bet$end[i]]
      gam <- c(alp, bet)
      
      rx <- as.matrix(int[, names(gam), drop = FALSE])
      
      tmp <- as.vector(exp(fx %*% the))
      yh <- tmp/(1+tmp)
      
      id <- c(id.alp$start[i]:id.alp$end[i], id.bet$start[i]:id.bet$end[i])
      gt[, id - offset] <- rx[, names(para)[id], drop = FALSE] * (fx0 * yh * (1-yh))
      
    }
    g.the[[j]] <- gt
    rm(gt)
  }
  
  g.the
  
}

