
# return optimal Sigma0, the covariance of auxiliary information
Sigma0.lm <- function(para, map, ref, model, nsample, outcome){
  
  #message('Estimating optimal covariance for auxiliary information...')
  
  ref$'(Intercept)' <- 1
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  n <- nrow(ref)
  
  lam <- para[map$lam]
  
  sigma <- para[map$the[1]]
  the <- para[map$the[-1]] # remove error term parameter
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  y <- ref[, outcome]
  
  hess <- matrix(0, nrow = nlam, ncol = nlam)
  score <- matrix(0, nrow = n, ncol = nlam)
  offset <- max(map$the)
  for(i in 1:nmodel){
    id.tau <- map$alp[[i]][1]
    tau <- para[id.tau]
    
    id.a <- alp.index.lm(map, i)
    alp.exist <- !is.null(id.a) # TRUE if at least one linear coef unknown; error term does not count
    if(alp.exist){
      alp <- para[id.a]
    }else{
      alp <- NULL
    }
    
    id.b <- map$bet[[i]]
    bet <- para[id.b]
    gam <- c(alp, bet)
    
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    if(alp.exist){
      ra <- as.matrix(rx[, names(alp), drop = FALSE])
    }else{
      ra <- NULL
    }
    rb <- as.matrix(rx[, names(bet), drop = FALSE])
    res <- as.vector(y - rx %*% gam)
    
    hess[id.tau - offset, id.tau - offset] <- 1/2/tau^2 - 1/tau^3 * mean(res^2)
    if(alp.exist){
      hess[id.tau - offset, id.a - offset] <- -1/tau^2 * t(t(ra) %*% res)/n
      hess[id.a - offset, id.tau - offset] <- t(hess[id.tau - offset, id.a - offset])
      hess[id.a - offset, id.a - offset] <- -1/tau * (t(ra) %*% ra)/n
    }
    
    hess[id.tau - offset, id.b - offset] <- -1/tau^2 * t(t(rb) %*% res)/n
    hess[id.b - offset, id.tau - offset] <- t(hess[id.tau - offset, id.b - offset])
    
    if(alp.exist){
      hess[id.a - offset, id.b - offset] <- -1/tau^2 * (t(ra) %*% rb)/n
      hess[id.b - offset, id.a - offset] <- t(hess[id.a - offset, id.b - offset])
    }
    
    hess[id.b - offset, id.b - offset] <- -1/tau * (t(rb) %*% rb)/n
    
    id <- c(id.tau, id.a, id.b)
    hess[id - offset, id - offset] <- nsample[i, i] * hess[id - offset, id - offset]
    
    ##############
    
    score[, id.tau - offset] <- -1/2/tau + 1/2/tau^2 * res^2
    if(alp.exist){
      score[, id.a - offset] <- 1/tau * (ra * res)
    }
    score[, id.b - offset] <- 1/tau * (rb * res)
    
    rm(id.tau, id.a, id.b, alp.exist)
    
  }
  
  info <- cov(score)
  for(i in 1:nmodel){
    for(j in i:nmodel){
      id1 <- c(map$alp[[i]], map$bet[[i]])
      id2 <- c(map$alp[[j]], map$bet[[j]])
      info[id1 - offset, id2 - offset] <- nsample[i, j] * info[id1 - offset, id2 - offset]
      info[id2 - offset, id1 - offset] <- t(info[id1 - offset, id2 - offset])
    }
  }
  
  V <- solve(hess) %*% info %*% solve(hess)
  id <- map$all.bet
  V <- V[id - offset, id - offset, drop = FALSE]
  colnames(V) <- names(para)[id]
  rownames(V) <- names(para)[id]
  
  V
  
}

