
# there is a problem with this function, why I assume outcome is available in ref??
# return optimal Sigma0, the covariance of auxiliary information
Sigma0.lo <- function(para, map, ref, model, nsample, outcome){
  
  #message('Estimating optimal covariance for auxiliary information...')
  
  ref$'(Intercept)' <- 1
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  n <- nrow(ref)
  
  lam <- para[map$lam]
  
  the <- para[map$the]
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  y <- ref[, outcome]
  
  hess <- matrix(0, nrow = nlam, ncol = nlam)
  score <- matrix(0, nrow = n, ncol = nlam)
  offset <- max(map$the)
  for(i in 1:nmodel){
    
    id.a <- alp.index.lo(map, i)
    alp.exist <- !is.null(id.a) # TRUE if at least one linear coef unknown (e.g. intercept)
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
    
    tmp <- as.vector(exp(rx %*% gam))
    yi <- tmp/(1+tmp)
    
    if(alp.exist){
      hess[id.a - offset, id.a - offset] <- - (t(ra) %*% (ra * yi * (1-yi)))/n
      hess[id.a - offset, id.b - offset] <- - (t(ra) %*% (rb * yi * (1-yi)))/n
      hess[id.b - offset, id.a - offset] <- t(hess[id.a - offset, id.b - offset])
    }
    
    hess[id.b - offset, id.b - offset] <- - (t(rb) %*% (rb * yi * (1-yi)))/n
    
    id <- c(id.a, id.b)
    hess[id - offset, id - offset] <- nsample[i, i] * hess[id - offset, id - offset]
    
    ##############
    
    if(alp.exist){
      score[, id.a - offset] <- ra * (y - yi)
    }
    score[, id.b - offset] <- rb * (y - yi)
    
    rm(id.a, id.b, alp.exist)
    
  }
  
  info <- cov(score)
  for(i in 1:nmodel){
    for(j in i:nmodel){
      # cannot use map$alp[[i]] as in Sigma0.lm
      # because for logistic model, map$alp[[i]] could be NA if even intercept is also provided
      # this could not happen in linear model because we assume that at least tau (error term) is not provided
      # alp.index.lo will convert NA to be NULL
      id1 <- c(alp.index.lo(map, i), map$bet[[i]])
      id2 <- c(alp.index.lo(map, j), map$bet[[j]])
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

