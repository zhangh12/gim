
# there is a problem with this function, why I assume outcome is available in ref??
# return optimal Sigma0, the covariance of auxiliary information
Sigma0.lm <- function(para, map, ref, model, nsample, tau, outcome){
  
  #message('Estimating optimal covariance for auxiliary information...')
  
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
    id.a <- alp.index.lm(map, i)
    alp.exist <- !is.null(id.a) # TRUE if at least one linear coef unknown
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
    
    if(alp.exist){
      hess[id.a - offset, id.a - offset] <- -1/tau[i] * (t(ra) %*% ra)/n
    }
    
    if(alp.exist){
      hess[id.a - offset, id.b - offset] <- -1/tau[i]^2 * (t(ra) %*% rb)/n
      hess[id.b - offset, id.a - offset] <- t(hess[id.a - offset, id.b - offset])
    }
    
    hess[id.b - offset, id.b - offset] <- -1/tau[i] * (t(rb) %*% rb)/n
    
    id <- c(id.a, id.b)
    hess[id - offset, id - offset] <- nsample[i, i] * hess[id - offset, id - offset]
    
    ##############
    
    if(alp.exist){
      score[, id.a - offset] <- 1/tau[i] * (ra * res)
    }
    score[, id.b - offset] <- 1/tau[i] * (rb * res)
    
    rm(id.a, id.b, alp.exist)
    
  }
  
  info <- cov(score)
  for(i in 1:nmodel){
    for(j in i:nmodel){
      id1 <- c(alp.index.lm(map, i), map$bet[[i]])
      id2 <- c(alp.index.lm(map, j), map$bet[[j]])
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

