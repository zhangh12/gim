

obj.lm <- function(para, para.id, map, data, ref, inv.V, bet0, outcome){
  
  data$'(Intercept)' <- 1
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  
  lam <- para[map$lam]
  
  sigma <- para[map$the[1]]
  the <- para[map$the[-1]]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  y <- data[, outcome]
  res <- y - as.vector(fx %*% the)
  
  g <- gfunction.lm(para, map, ref)
  
  pr <- as.vector(1/(1+g %*% lam))
  
  np <- length(para)
  n <- nrow(data)
  
  bet <- para[map$all.bet]
  
  plogL <- -n/2 * log(sigma) - 1/2/sigma * sum(res^2) + sum(log(pr)) -1/2 * as.vector(t(bet0 - bet) %*% inv.V %*% (bet0 - bet))
  names(plogL) <- NULL
  plogL
  
}

