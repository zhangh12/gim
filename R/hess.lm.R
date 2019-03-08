
hess.lm <- function(para, map, data, ref, inv.V, bet0, outcome){
  
  if(0){
    h <- numDeriv::jacobian(score.lm, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome)
    colnames(h) <- names(para)
    rownames(h) <- names(para)
    return(h)
  }
  
  np <- length(para)
  h <- matrix(0, np, np)
  
  nmodel <- length(map$bet)
  nlam <- length(map$lam)
  n <- nrow(data)
  
  lam <- para[map$lam]
  tau <- para[map$all.tau]
  
  sigma <- para[map$the[1]]
  the <- para[map$the[-1]]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  y <- data[, outcome]
  res <- y - as.vector(fx %*% the)
  
  g <- gfunction.lm(para, map, ref)
  g.the <- gfunction.the.lm(para, map, ref)
  g.alp <- gfunction.alp.lm(para, map, ref)
  g.bet <- gfunction.bet.lm(para, map, ref)
  
  g.the.lam <- gfunction.the.lam.lm(g.the, lam)
  g.alp.lam <- gfunction.alp.lam.lm(g.alp, lam)
  g.bet.lam <- gfunction.bet.lam.lm(g.bet, lam)
  
  pr <- as.vector(1/(1+g %*% lam))
  pr2 <- pr^2
  
  ########## ell.lam.lam, done
  
  h[map$lam, map$lam] <- t(g) %*% (g * pr2)
  
  ########## ell.lam.the (first one is for sigma), done
  
  nthe <- length(map$the)
  for(i in 1:nthe){
    h[map$lam, map$the[i]] <- - t(g.the[[i]]) %*% pr + t(g) %*% (g.the.lam[, i] * pr2)
  }
  h[map$the, map$lam] <- t(h[map$lam, map$the])
  
  
  if(!is.null(map$all.alp)){
    ########## ell.lam.alp, done
    nalp <- length(map$all.alp)
    offset <- max(map$the)
    for(i in 1:nalp){
      h[map$lam, map$all.alp[i]] <- -t(g.alp[[i]]) %*% pr + t(g) %*% (g.alp.lam[, i] * pr2)
    }
    h[map$all.alp, map$lam] <- t(h[map$lam, map$all.alp])
    
    ########## ell.the.alp (inc ell.sigma.alp), done
    nthe <- length(map$the)
    for(i in 1:nthe){
      k <- 0
      for(j in map$all.alp){
        k <- k + 1
        h[map$the[i], j] <- sum(g.the.lam[, i] * g.alp.lam[, k] * pr2)
        h[j, map$the[i]] <- h[map$the[i], j]
      }
    }
  }
  
  ########## ell.lam.bet, done
  
  nbet <- length(map$all.bet)
  offset <- max(map$the)
  for(i in 1:nbet){
    h[map$lam, map$all.bet[i]] <- -t(g.bet[[i]]) %*% pr + t(g) %*% (g.bet.lam[, i] * pr2)
  }
  h[map$all.bet, map$lam] <- t(h[map$lam, map$all.bet])
  
  ########## ell.the.the (not include sigma), done
  # ell.sigma.sigma
  h[map$the[1], map$the[1]] <- n/2/sigma^2 - 1/sigma^3 * sum(res^2)
  
  # ell.sigma.the
  h[map$the[1], map$the[-1]] <- -1/sigma^2 * t(t(fx) %*% res)
  h[map$the[-1], map$the[1]] <- t(h[map$the[1], map$the[-1]])
  
  # ell.the.the
  h[map$the[-1], map$the[-1]] <- -1/sigma * (t(fx) %*% fx)
  nthe <- length(map$the)
  for(j in 1:nthe){
    id.j <- map$the[j]
    for(l in j:nthe){
      id.l <- map$the[l]
      h[id.j, id.l] <- h[id.j, id.l] + sum(g.the.lam[, j] * g.the.lam[, l] * pr2)
      h[id.l, id.j] <- t(h[id.j, id.l])
    }
  }
  
  ########## ell.the.bet, include ell.sigma.bet, done
  nthe <- length(map$the)
  for(i in 1:nthe){
    k <- 0
    for(j in map$all.bet){
      k <- k + 1
      h[map$the[i], j] <- sum(g.the.lam[, i] * g.bet.lam[, k] * pr2)
      h[j, map$the[i]] <- h[map$the[i], j]
    }
  }
  
  ########## ell.alp.alp, done
  id.map.alp <- list()
  k <- 0
  for(i in 1:nmodel){
    id.a <- alp.index.lm(map, i)
    if(is.null(id.a)){
      next
    }
    
    for(j in id.a){
      k <- k + 1
      id.map.alp[[j]] <- k
    }
  }
  
  for(j in map$all.alp){
    for(l in map$all.alp){
      if(l < j){
        next
      }
      
      h[j, l] <- sum(g.alp.lam[, id.map.alp[[j]]] * g.alp.lam[, id.map.alp[[l]]] * pr2)
      h[l, j] <- h[j, l]
    }
  }
  
  ########## ell.alp.bet, done
  k <- 0
  for(j in map$all.bet){
    k <- k + 1
    for(l in map$all.alp){
      h[j, l] <- sum(g.alp.lam[, id.map.alp[[l]]] * g.bet.lam[, k] * pr2)
      h[l, j] <- h[j, l]
    }
  }
  
  ########## ell.bet.bet
  id.map.bet <- list()
  k <- 0
  for(i in 1:nmodel){
    id.b <- map$bet[[i]]
    for(j in id.b){
      k <- k + 1
      id.map.bet[[j]] <- k
    }
  }
  
  for(j in map$all.bet){
    for(l in map$all.bet){
      if(l < j){
        next
      }
      h[j, l] <- sum(g.bet.lam[, id.map.bet[[j]]] * g.bet.lam[, id.map.bet[[l]]] * pr2)
      h[l, j] <- h[j, l]
    }
  }
  
  h[map$all.bet, map$all.bet] <- h[map$all.bet, map$all.bet] - inv.V
  
  
  ##########
  
  colnames(h) <- names(para)
  rownames(h) <- names(para)
  
  h
  
}
