
# return optimal Sigma0, the covariance of auxiliary information
Sigma0.lm <- function(para, para.id, int, model, nsample, outcome = 'y'){
  
  message('Estimating optimal covariance for auxiliary information...')
  
  int$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.alp)
  nt <- max(id.lam)
  n <- nrow(int)
  
  lam <- para[id.lam$start[1]:id.lam$end[1]]
  
  sigma <- para[id.the$start[1]]
  the <- para[(id.the$start[1]+1):id.the$end[1]]
  fx <- as.matrix(int[, names(the), drop = FALSE])
  y <- int[, outcome]
  
  hess <- matrix(0, nrow = nt, ncol = nt)
  score <- matrix(0, nrow = n, ncol = nt)
  offset <- max(id.the)
  for(i in 1:nmodel){
    tau <- para[id.alp$start[i]]
    alp <- para[(id.alp$start[i] + 1):id.alp$end[i]]
    bet <- para[id.bet$start[i]:id.bet$end[i]]
    gam <- c(alp, bet)
    
    rx <- as.matrix(int[, names(gam), drop = FALSE])
    ra <- as.matrix(rx[, names(alp), drop = FALSE])
    rb <- as.matrix(rx[, names(bet), drop = FALSE])
    res <- as.vector(y - rx %*% gam)
    
    hess[id.alp$start[i] - offset, id.alp$start[i] - offset] <- 1/2/tau^2 - 1/tau^3 * mean(res^2)
    id.a <- (id.alp$start[i]+1):id.alp$end[i]
    hess[id.alp$start[i] - offset, id.a - offset] <- -1/tau^2 * t(t(ra) %*% res)/n
    hess[id.a - offset, id.alp$start[i] - offset] <- t(hess[id.alp$start[i] - offset, id.a - offset])
    
    hess[id.a - offset, id.a - offset] <- -1/tau * (t(ra) %*% ra)/n
    
    id.b <- id.bet$start[i]:id.bet$end[i]
    hess[id.alp$start[i] - offset, id.b - offset] <- -1/tau^2 * t(t(rb) %*% res)/n
    hess[id.b - offset, id.alp$start[i] - offset] <- t(hess[id.alp$start[i] - offset, id.b - offset])
    
    hess[id.a - offset, id.b - offset] <- -1/tau^2 * (t(ra) %*% rb)/n
    hess[id.b - offset, id.a - offset] <- t(hess[id.a - offset, id.b - offset])
    
    hess[id.b - offset, id.b - offset] <- -1/tau * (t(rb) %*% rb)/n
    
    id <- c(id.alp$start[i], id.a, id.b)
    hess[id - offset, id - offset] <- nsample[i, i] * hess[id - offset, id - offset]
    
    ##############
    
    score[, id.alp$start[i] - offset] <- -1/2/tau + 1/2/tau^2 * res^2
    score[, id.a - offset] <- 1/tau * (ra * res)
    score[, id.b - offset] <- 1/tau * (rb * res)
    
  }
  
  info <- cov(score)
  for(i in 1:nmodel){
    for(j in i:nmodel){
      id1 <- c(id.alp$start[i]:id.alp$end[i], id.bet$start[i]:id.bet$end[i])
      id2 <- c(id.alp$start[j]:id.alp$end[j], id.bet$start[j]:id.bet$end[j])
      info[id1 - offset, id2 - offset] <- nsample[i, j] * info[id1 - offset, id2 - offset]
      info[id2 - offset, id1 - offset] <- t(info[id1 - offset, id2 - offset])
    }
  }
  
  V <- solve(hess) %*% info %*% solve(hess)
  id <- min(id.bet):max(id.bet) - offset
  V <- V[id, id, drop = FALSE]
  
  V
  
}

