

obj.lo <- function(para, map, data, ref, inv.V, bet0, outcome){
  
  data$'(Intercept)' <- 1
  
  nmodel <- length(map$bet)
  
  lam <- para[map$lam]
  the <- para[map$the]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  y <- data[, outcome]
  lin <- as.vector(fx %*% the)
  
  g <- gfunction.lo(para, map, ref)
  
  pr <- as.vector(1/(1+g %*% lam))
  
  bet <- para[map$all.bet]
  
  plogL <- sum(y*lin) - sum(log(1+exp(lin))) + sum(log(pr)) -1/2 * as.vector(t(bet0 - bet) %*% inv.V %*% (bet0 - bet))
  names(plogL) <- NULL
  plogL
  
}

