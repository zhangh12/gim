

obj.lo <- function(para, para.id, int, inv.V, bet0, outcome = 'y'){
  
  int$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.alp)
  nt <- max(id.lam)
  
  lam <- para[id.lam$start[1]:id.lam$end[1]]
  
  the <- para[id.the$start[1]:id.the$end[1]]
  fx <- as.matrix(int[, names(the), drop = FALSE])
  y <- int[, outcome]
  lin <- as.vector(fx %*% the)
  
  g <- gfunction.lo(para, para.id, int)
  
  pr <- as.vector(1/(1+g %*% lam))
  
  np <- length(para)
  n <- nrow(int)
  
  bet <- para[min(id.bet):max(id.bet)]
  
  plogL <- sum(y*lin) - sum(log(1+exp(lin))) + sum(log(pr)) -1/2 * as.vector(t(bet0 - bet) %*% inv.V %*% (bet0 - bet))
  names(plogL) <- NULL
  plogL
  
}

