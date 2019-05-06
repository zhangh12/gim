
# return optimal Sigma0, the covariance of auxiliary information
Sigma0.lm <- function(para, map, ref, model, nsample, outcome){
  
  #message('Estimating optimal covariance for auxiliary information...')
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  n <- nrow(ref)
  
  sigma <- para[map$the[1]]
  the <- para[map$the[-1]]
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  x.the <- as.vector(fx %*% the)
  
  hess <- matrix(0, nrow = nlam, ncol = nlam)
  info <- matrix(0, nrow = nlam, ncol = nlam)
  offset <- max(map$the)
  
  for(i in 1:nmodel){
    id1 <- c(alp.index.lm(map, i), map$bet[[i]])
    gam1 <- para[id1]
    rx1 <- as.matrix(ref[, names(gam1), drop = FALSE])
    hess[id1 - offset, id1 - offset] <- - nsample[i, i] * (t(rx1) %*% rx1) / n
    x.gam1 <- as.vector(rx1 %*% gam1)
    r1 <- x.the - x.gam1
    for(j in i:nmodel){
      id2 <- c(alp.index.lm(map, j), map$bet[[j]])
      gam2 <- para[id2]
      rx2 <- as.matrix(ref[, names(gam2), drop = FALSE])
      x.gam2 <- as.vector(rx2 %*% gam2)
      r2 <- x.the - x.gam2
      tmp <- sigma * (t(rx1) %*% rx2) / n + t(rx1) %*% (rx2 * (r1 * r2)) / n
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

