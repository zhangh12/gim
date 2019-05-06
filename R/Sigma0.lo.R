
# there is a problem with this function, why I assume outcome is available in ref??
# return optimal Sigma0, the covariance of auxiliary information
Sigma0.lo <- function(para, map, ref, model, nsample, outcome){
  
  #message('Estimating optimal covariance for auxiliary information...')
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  n <- nrow(ref)
  
  the <- para[map$the]
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  exp.x.the <- exp(as.vector(fx %*% the))
  y <- exp.x.the / (1 + exp.x.the)
  
  hess <- matrix(0, nrow = nlam, ncol = nlam)
  info <- matrix(0, nrow = nlam, ncol = nlam)
  offset <- max(map$the)
  for(i in 1:nmodel){
    id1 <- c(alp.index.lo(map, i), map$bet[[i]])
    gam1 <- para[id1]
    rx1 <- as.matrix(ref[, names(gam1), drop = FALSE])
    exp.x.gam1 <- exp(as.vector(rx1 %*% gam1))
    y1 <- exp.x.gam1 / (1 + exp.x.gam1)
    hess[id1 - offset, id1 - offset] <- -nsample[i, i] * (t(rx1) %*% (rx1 * y1 * (1 - y1))) / n
    for(j in i:nmodel){
      id2 <- c(alp.index.lo(map, j), map$bet[[j]])
      gam2 <- para[id2]
      rx2 <- as.matrix(ref[, names(gam2), drop = FALSE])
      exp.x.gam2 <- exp(as.vector(rx2 %*% gam2))
      y2 <- exp.x.gam2 / (1 + exp.x.gam2)
      tmp <- t(rx1) %*% (rx2 * (y * (1 - y1 - y2) + y1 * y2)) / n
      info[id1 - offset, id2 - offset] <- nsample[i, j] * tmp
      info[id2 - offset, id1 - offset] <- t(info[id1 - offset, id2 - offset])
      rm(rx2)
    }
    rm(rx1)
  }
  
  V <- solve(hess) %*% info %*% solve(hess)
  id <- map$all.bet
  V <- V[id - offset, id - offset, drop = FALSE]
  colnames(V) <- names(para)[id]
  rownames(V) <- names(para)[id]
  
  V
  
}

