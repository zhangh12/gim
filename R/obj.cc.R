

obj.cc <- function(para, map, data, ref, inv.V, bet0, sample.info, outcome){
  
  ncase <- sample.info$ncase
  nctrl <- sample.info$nctrl
  
  nmodel <- length(map$bet)
  lam <- para[map$lam]
  the <- para[map$the]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  y <- data[, outcome]
  
  tilt <- tilt.cc(para, map, ref)
  Delta <- tilt$Delta
  delta <- tilt$delta
  
  g <- gfunction.cc(para, map, ref, Delta, delta, ncase, nctrl)
  
  pr <- as.vector(1/(1+g %*% lam))
  bet <- para[map$all.bet]
  
  plogL <- sum(y * log(Delta)) + sum(log(pr)) -1/2 * as.vector(t(bet0 - bet) %*% inv.V %*% (bet0 - bet))
  names(plogL) <- NULL
  plogL
  
}

