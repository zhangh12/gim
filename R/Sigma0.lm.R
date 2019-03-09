
# return optimal Sigma0, the covariance of auxiliary information
Sigma0.lm <- function(para, map, ref, model, nsample, outcome){
  
  #message('Estimating optimal covariance for auxiliary information...')
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  n <- nrow(ref)
  
  sigma <- para[map$the[1]]
  
  hess <- matrix(0, nrow = nlam, ncol = nlam)
  info <- matrix(0, nrow = nlam, ncol = nlam)
  offset <- max(map$the)
  
  for(i in 1:nmodel){
    id1 <- c(alp.index.lm(map, i), map$bet[[i]])
    gam1 <- para[id1]
    rx1 <- as.matrix(ref[, names(gam1), drop = FALSE])
    hess[id1 - offset, id1 - offset] <- - nsample[i, i] * (t(rx1) %*% rx1) / n
    for(j in i:nmodel){
      id2 <- c(alp.index.lm(map, j), map$bet[[j]])
      gam2 <- para[id2]
      rx2 <- as.matrix(ref[, names(gam2), drop = FALSE])
      info[id1 - offset, id2 - offset] <- nsample[i, j] * (t(rx1) %*% rx2) / n
      info[id2 - offset, id1 - offset] <- t(info[id1 - offset, id2 - offset])
      rm(rx2)
    }
    rm(rx1)
  }
  
  info <- sigma * info
  
  V <- solve(hess) %*% info %*% solve(hess)
  id <- map$all.bet
  V <- V[id - offset, id - offset, drop = FALSE]
  colnames(V) <- names(para)[id]
  rownames(V) <- names(para)[id]
  
  V
  
}

