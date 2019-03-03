
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
  
  ########## ell.lam.lam
  
  h[map$lam, map$lam] <- t(g) %*% (g * pr2)
  
  ########## ell.lam.the (first one is for sigma)
  
  nthe <- length(map$the)
  for(i in 1:nthe){
    h[map$lam, map$the[i]] <- -t(g.the[[i]]) %*% pr + t(g) %*% (g.the.lam[, i] * pr2)
  }
  h[map$the, map$lam] <- t(h[map$lam, map$the])
  
  ########## ell.lam.alp
  
  # nalp > 0 as tau is in alp
  nalp <- length(map$all.alp)
  offset <- max(map$the)
  for(i in 1:nalp){
    h[map$lam, map$all.alp[i]] <- -t(g.alp[[i]]) %*% pr + t(g) %*% (g.alp.lam[, i] * pr2)
  }
  h[map$all.alp, map$lam] <- t(h[map$lam, map$all.alp])
  
  ########## ell.lam.bet
  
  nbet <- length(map$all.bet)
  offset <- max(map$the)
  for(i in 1:nbet){
    h[map$lam, map$all.bet[i]] <- -t(g.bet[[i]]) %*% pr + t(g) %*% (g.bet.lam[, i] * pr2)
    h[map$all.bet[i], map$lam] <- t(h[map$lam, map$all.bet[i]])
  }
  h[map$all.bet, map$lam] <- t(h[map$lam, map$all.bet])
  
  ########## ell.sigma.sigma
  
  h[map$the[1], map$the[1]] <- n/2/sigma^2 - 1/sigma^3 * sum(res^2) + sum(g.the.lam[, 1]^2 * pr2)
  
  ########## ell.sigma.the
  
  nthe <- length(map$the)
  h[map$the[1], map$the[-1]] <- -1/sigma^2 * t(t(fx) %*% res)
  for(i in 2:nthe){
    h[map$the[1], map$the[i]] <- h[map$the[1], map$the[i]] + sum(g.the.lam[, 1] * g.the.lam[, i] * pr2)
    h[map$the[i], map$the[1]] <- t(h[map$the[1], map$the[i]])
  }
  
  ########## ell.sigma.alp (first one is for tau)
  
  offset <- max(map$the)
  for(j in map$all.alp){
    h[map$the[1], j] <- sum(g.the.lam[, 1] * g.alp.lam[, j - offset] * pr2)
    h[j, map$the[1]] <- h[map$the[1], j]
  }
  
  ########## ell.sigma.bet
  
  offset <- max(map$all.alp)
  for(j in map$all.bet){
    h[map$the[1], j] <- sum(g.the.lam[, 1] * g.bet.lam[, j - offset] * pr2)
    h[j, map$the[1]] <- h[map$the[1], j]
  }
  
  ########## ell.the.the (not include sigma)
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  g.the2 <- gfunction.the2.lm(para, map, ref, lam)
  
  if(0){
    foo1 <- function(the.l, para, map, ref, lam, j, l){
      
      para[map$the[l]] <- the.l
      sum(gfunction.the.lm(para, map, ref)[[j]] %*% lam)
      
    }
  }
  
  if(0){
    for(j in 2:length(map$the)){
      for(l in j:length(map$the)){
        tt <- grad(foo1, para[map$the[l]], para=para, map=map, ref=ref, lam=lam, j = j, l = l)
        cat(j,'\t', l, '\t', tt, '\t', sum(g.the2[[foo(j,l)]]), '\t', tt-sum(g.the2[[foo(j,l)]]), '\n')
      }
    }
  }
  
  
  h[map$the[-1], map$the[-1]] <- -1/sigma * (t(fx) %*% fx)
  nthe <- length(map$the)
  for(j in 2:nthe){
    id.j <- map$the[j]
    for(l in j:nthe){
      id.l <- map$the[l]
      h[id.j, id.l] <- h[id.j, id.l] - 
        t(g.the2[[foo(j,l)]]) %*% pr + sum(g.the.lam[, j] * g.the.lam[, l] * pr2)
      h[id.l, id.j] <- t(h[id.j, id.l])
    }
  }
  
  ########## ell.the.alp (include tau)
  
  g.the.alp <- gfunction.the.alp.lm(para, map, ref, lam)
  
  nthe <- length(map$the)
  
  for(j in 2:nthe){
    id.j <- map$the[j]
    k <- 0
    for(i in 1:nmodel){
      id.tau <- map$alp[[i]][1]
      k <- k + 1
      h[id.j, id.tau] <- sum(g.the.lam[, j] * g.alp.lam[, k] * pr2)
      h[id.tau, id.j] <- h[id.j, id.tau]
      id.a <- alp.index.lm(map, i)
      if(is.null(id.a)){
        next
      }
      
      for(l in id.a){
        k <- k + 1
        h[id.j, l] <- -sum(g.the.alp[[foo(j,l)]] * pr) + sum(g.the.lam[, j] * g.alp.lam[, k] * pr2)
        h[l, id.j] <- t(h[id.j, l])
      }
    }
  }
  
  ########## ell.the.bet
  
  g.the.bet <- gfunction.the.bet.lm(para, map, ref, lam)
  
  nthe <- length(map$the)
  
  for(j in 2:nthe){
    id.j <- map$the[j]
    k <- 0
    for(i in 1:nmodel){
      id.tau <- map$alp[[i]][1]
      id.b <- map$bet[[i]]
      
      for(l in id.b){
        k <- k + 1
        h[id.j, l] <- -sum(g.the.bet[[foo(j,l)]] * pr) + sum(g.the.lam[, j] * g.bet.lam[, k] * pr2)
        h[l, id.j] <- t(h[id.j, l])
      }
    }
  }
  
  ########## ell.alp.alp
  
  g.alp2 <- gfunction.alp2.lm(para, map, ref, lam)
  
  offset <- max(map$the)
  id.map.alp <- list()
  k <- 0
  for(i in 1:nmodel){
    id.tau <- map$alp[[i]][1]
    k <- k + 1
    id.map.alp[[id.tau]] <- k
    id.a <- alp.index.lm(map, i)
    if(is.null(id.a)){
      next
    }
    
    for(j in id.a){
      k <- k + 1
      id.map.alp[[j]] <- k
    }
  }
  
  # ell.tau.tau
  for(i1 in 1:nmodel){
    id.tau1 <- map$alp[[i1]][1]
    for(i2 in i1:nmodel){
      id.tau2 <- map$alp[[i2]][1]
      
      h[id.tau1, id.tau2] <- sum(g.alp.lam[, id.map.alp[[id.tau1]]] * 
                                   g.alp.lam[, id.map.alp[[id.tau2]]] * pr2)
      h[id.tau2, id.tau1] <- h[id.tau1, id.tau2]
    }
  }
  
  id.alp <- setdiff(map$all.alp, map$all.tau)
  
  # ell.alp.tau and ell.alp.alp
  for(j in id.alp){
    for(l in map$all.tau){
      h[l, j] <- sum(g.alp.lam[, id.map.alp[[l]]] * g.alp.lam[, id.map.alp[[j]]] * pr2)
      h[j, l] <- h[l, j]
    }
    
    for(l in id.alp){
      if(l < j){
        next
      }
      
      h[j, l] <- -sum(g.alp2[[foo(j,l)]] * pr) + 
        sum(g.alp.lam[, id.map.alp[[j]]] * g.alp.lam[, id.map.alp[[l]]] * pr2)
      h[l, j] <- h[j, l]
    }
  }
  
  ########## ell.tau.bet
  
  for(j in map$all.tau){
    k <- 0
    for(l in map$all.bet){
      k <- k + 1
      h[j, l] <- sum(g.alp.lam[, id.map.alp[[j]]] * g.bet.lam[, k] * pr2)
      h[l, j] <- h[j, l]
    }
  }
  
  ########## ell.alp.bet
  
  g.alp.bet <- gfunction.alp.bet.lm(para, map, ref, lam)
  
  offset <- max(map$the)
  id.alp <- setdiff(map$all.alp, map$all.tau)
  k <- 0
  for(j in map$all.bet){
    k <- k + 1
    for(l in id.alp){
      h[j, l] <- -sum(g.alp.bet[[foo(j,l)]] * pr) + 
        sum(g.alp.lam[, id.map.alp[[l]]] * g.bet.lam[, k] * pr2)
      h[l, j] <- h[j, l]
    }
  }
  
  ########## ell.bet.bet
  
  g.bet2 <- gfunction.bet2.lm(para, map, ref, lam)
  id.map.bet <- list()
  k <- 0
  for(i in 1:nmodel){
    id.b <- map$bet[[i]]
    for(j in id.b){
      k <- k + 1
      id.map.bet[[j]] <- k
      for(l in id.b){
        if(l < j){
          next
        }
        
        h[j, l] <- -sum(g.bet2[[foo(j,l)]] * pr)
        h[l, j] <- h[j, l]
      }
    }
  }
  
  for(j in map$all.bet){
    for(l in map$all.bet){
      if(l < j){
        next
      }
      h[j, l] <- h[j, l] + sum(g.bet.lam[, id.map.bet[[j]]] * g.bet.lam[, id.map.bet[[l]]] * pr2)
      h[l, j] <- h[j, l]
    }
  }
  
  h[map$all.bet, map$all.bet] <- h[map$all.bet, map$all.bet] - inv.V
  
  
  ##########
  
  colnames(h) <- names(para)
  rownames(h) <- names(para)
  h
  
}
