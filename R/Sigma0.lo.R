
# return optimal Sigma0, the covariance of auxiliary information
Sigma0.lo <- function(para, para.id, data, model, nsample, outcome = 'y'){
  
  #message('Estimating optimal covariance for auxiliary information...')
  
  int$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.bet)
  nlam <- max(id.lam)
  n <- nrow(int)
  
  lam <- para[id.lam$start[1]:id.lam$end[1]]
  
  the <- para[id.the$start[1]:id.the$end[1]]
  fx <- as.matrix(int[, names(the), drop = FALSE])
  y <- int[, outcome]
  
  hess <- matrix(0, nrow = nlam, ncol = nlam)
  score <- matrix(0, nrow = n, ncol = nlam)
  offset <- max(id.the)
  for(i in 1:nmodel){
    
    id.a <- alp.index.lo(id.alp, i)
    alp.exist <- !is.null(id.a)
    if(alp.exist){
      alp <- para[id.a]
    }else{
      alp <- NULL
    }
    
    id.b <- id.bet$start[i]:id.bet$end[i]
    bet <- para[id.b]
    gam <- c(alp, bet)
    
    rx <- as.matrix(int[, names(gam), drop = FALSE])
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
      id1 <- c(alp.index.lo(id.alp, i), id.bet$start[i]:id.bet$end[i])
      id2 <- c(alp.index.lo(id.alp, j), id.bet$start[j]:id.bet$end[j])
      info[id1 - offset, id2 - offset] <- nsample[i, j] * info[id1 - offset, id2 - offset]
      info[id2 - offset, id1 - offset] <- t(info[id1 - offset, id2 - offset])
    }
  }
  
  V <- solve(hess) %*% info %*% solve(hess)
  id <- min(id.bet):max(id.bet)
  V <- V[id - offset, id - offset, drop = FALSE]
  colnames(V) <- names(para)[id]
  rownames(V) <- names(para)[id]
  
  V
  
}

