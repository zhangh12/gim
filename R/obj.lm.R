

obj.lm <- function(para, para.id, data, inv.V, bet0, outcome = 'y'){
  
  data$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.bet)
  nlam <- max(id.lam)
  
  lam <- para[id.lam$start[1]:id.lam$end[1]]
  
  sigma <- para[id.the$start[1]]
  the <- para[(id.the$start[1]+1):id.the$end[1]]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  y <- data[, outcome]
  res <- y - as.vector(fx %*% the)
  
  g <- gfunction.lm(para, para.id, data)
  
  pr <- as.vector(1/(1+g %*% lam))
  
  np <- length(para)
  n <- nrow(data)
  
  bet <- para[min(id.bet):max(id.bet)]
  
  plogL <- -n/2 * log(sigma) - 1/2/sigma * sum(res^2) + sum(log(pr)) -1/2 * as.vector(t(bet0 - bet) %*% inv.V %*% (bet0 - bet))
  names(plogL) <- NULL
  plogL
  
}

