
# return optimal Sigma0, the covariance of auxiliary information
Sigma0.cc <- function(para, map, ref, model, sample.info, pr0, Delta, outcome){
  
  #message('Estimating optimal covariance for auxiliary information...')
  
  ncase <- sample.info$ncase
  nctrl <- sample.info$nctrl
  
  ud <- update.pr0.Delta(para, map, ref, Delta, ncase, nctrl, pr0)
  pr0 <- ud$pr0
  Delta <- ud$Delta
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  
  hess <- matrix(0, nrow = nlam - 1, ncol = nlam - 1)
  info <- hess
  offset <- max(map$the)
  
  E0 <- rep(list(NULL), nmodel)
  for(i in 1:nmodel){
    
    id.i <- c(alp.index.cc(map, i), map$bet[[i]])
    gam <- para[id.i]
    
    rx.i <- as.matrix(ref[, names(gam), drop = FALSE])
    
    delta.i <- as.vector(exp(rx.i %*% gam))
    n1.i <- diag(ncase)[i]
    n0.i <- diag(nctrl)[i]
    rho.i <- n1.i/n0.i
    tmp <- -n1.i * delta.i * (1 + rho.i * Delta) / (1 + rho.i * delta.i)^2 * pr0
    
    hess[id.i - offset, id.i - offset] <- (t(rx.i) %*% (rx.i * tmp))
    
    for(j in i:nmodel){
      
      id.j <- c(alp.index.lo(map, j), map$bet[[j]])
      gam <- para[id.j]
      
      rx.j <- as.matrix(ref[, names(gam), drop = FALSE])
      
      delta.j <- as.vector(exp(rx.j %*% gam))
      n1.j <- diag(ncase)[j]
      n0.j <- diag(nctrl)[j]
      rho.j <- n1.j/n0.j
      
      n1.ij <- ncase[i, j]
      n0.ij <- nctrl[i, j]
      tmp <- (n1.ij * Delta + n0.ij * rho.i * rho.j * delta.i * delta.j) / (1 + rho.i * delta.i) / (1 + rho.j * delta.j) * pr0
      
      if(is.null(E0[[i]])){
        E0[[i]] <- t(rx.i) %*% (Delta / (1 + rho.i * delta.i) * pr0)
      }
      
      if(is.null(E0[[j]])){
        E0[[j]] <- t(rx.j) %*% (Delta / (1 + rho.j * delta.j) * pr0)
      }
      
      info[id.i - offset, id.j - offset] <- (t(rx.i) %*% (rx.j * tmp))
      info[id.i - offset, id.j - offset] <- info[id.i - offset, id.j - offset] - (n1.ij + n0.ij * rho.i * rho.j) * (E0[[i]] %*% t(E0[[j]]))
      
      if(i == j){
        next
      }
      
      info[id.j - offset, id.i - offset] <- t(info[id.i - offset, id.j - offset])
    }
    
  }
  
  V <- solve(hess) %*% info %*% solve(hess)
  id <- map$all.bet
  V <- V[id - offset, id - offset, drop = FALSE]
  colnames(V) <- names(para)[id]
  rownames(V) <- names(para)[id]
  
  V
  
}

