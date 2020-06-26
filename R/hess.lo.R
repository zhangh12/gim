
hess.lo <- function(para, map, data, ref, inv.V, bet0, outcome){
  
  if(0){
    h <- numDeriv::jacobian(score.lo, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome)
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
  the <- para[map$the]
  
  fx <- as.matrix(data[, names(the), drop = FALSE])
  tilt <- tilt.cc(para, map, ref)
  Delta <- tilt$Delta
  delta <- tilt$delta
  yh <- Delta/(1 + Delta)
  
  g <- gfunction.lo(para, map, ref)
  g.the <- gfunction.the.lo(para, map, ref)
  g.alp <- gfunction.alp.lo(para, map, ref)
  g.bet <- gfunction.bet.lo(para, map, ref)
  
  g.the.lam <- gfunction.the.lam.lo(g.the, lam)
  g.alp.lam <- gfunction.alp.lam.lo(g.alp, lam)
  g.bet.lam <- gfunction.bet.lam.lo(g.bet, lam)
  
  pr <- as.vector(1/(1+g %*% lam))
  pr2 <- pr^2
  
  ########## ell.lam.lam
  
  h[map$lam, map$lam] <- t(g) %*% (g * pr2)
  
  ########## ell.lam.the
  
  nthe <- length(map$the)
  for(i in 1:nthe){
    h[map$lam, map$the[i]] <- - t(g.the[[i]]) %*% pr + t(g) %*% (g.the.lam[, i] * pr2)
  }
  h[map$the, map$lam] <- t(h[map$lam, map$the])
  
  if(!is.null(map$all.alp)){
    ########## ell.lam.alp
    nalp <- length(map$all.alp)
    for(i in 1:nalp){
      h[map$lam, map$all.alp[i]] <- -t(g.alp[[i]]) %*% pr + t(g) %*% (g.alp.lam[, i] * pr2)
    }
    h[map$all.alp, map$lam] <- t(h[map$lam, map$all.alp])
    
    ########## ell.the.alp
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
  
  ########## ell.lam.bet
  
  nbet <- length(map$all.bet)
  offset <- max(map$the)
  for(i in 1:nbet){
    h[map$lam, map$all.bet[i]] <- -t(g.bet[[i]]) %*% pr + t(g) %*% (g.bet.lam[, i] * pr2)
  }
  h[map$all.bet, map$lam] <- t(h[map$lam, map$all.bet])
  
  ########## ell.the.the
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  g.the2 <- gfunction.the2.lo(para, map, ref, Delta, lam, pr)
  h[map$the, map$the] <- -t(fx) %*% (fx * yh * (1 - yh))
  nthe <- length(map$the)
  for(j in 1:nthe){
    id.j <- map$the[j]
    for(l in j:nthe){
      id.l <- map$the[l]
      h[id.j, id.l] <- h[id.j, id.l] - g.the2[[foo(j,l)]] + 
        sum(g.the.lam[, j] * g.the.lam[, l] * pr2)
      h[id.l, id.j] <- t(h[id.j, id.l])
    }
  }
  
  ########## ell.the.bet
  nthe <- length(map$the)
  for(i in 1:nthe){
    k <- 0
    for(j in map$all.bet){
      k <- k + 1
      h[map$the[i], j] <- sum(g.the.lam[, i] * g.bet.lam[, k] * pr2)
      h[j, map$the[i]] <- h[map$the[i], j]
    }
  }
  
  ####### ell.alp.alp
  if(!is.null(map$all.alp)){
    g.alp2 <- gfunction.alp2.lo(para, map, ref, delta, lam, pr)
    
    kk <- seq_along(map$all.alp)
    names(kk) <- as.character(map$all.alp)
    k <- 0
    for(i in 1:nmodel){
      id.a <- alp.index.lo(map, i)
      if(is.null(id.a)){
        next
      }
      
      for(j in id.a){
        k <- k + 1
        tmp <- g.alp.lam[, k]
        for(l in id.a){
          if(l < j){
            next
          }
          
          h[j, l] <- -g.alp2[[foo(j, l)]] + 
            sum(tmp * g.alp.lam[, kk[as.character(l)]] * pr2)
          h[l, j] <- t(h[j, l])
        }
      }
    }
    
    rm(g.alp2)
    gc()
    
    k <- 0
    for(j in map$all.alp){
      k <- k + 1
      tmp <- g.alp.lam[, k]
      for(l in map$all.alp){
        if(l <= j){
          next
        }
        
        if(!is.na(h[j, l])){
          next
        }
        
        h[j, l] <- sum(tmp * g.alp.lam[, kk[as.character(l)]] * pr2)
        h[l, j] <- t(h[j, l])
      }
    }
  }
  ####### ell.alp.bet
  
  if(!is.null(map$all.alp)){
    g.alp.bet <- gfunction.alp.bet.lo(para, map, ref, delta, lam, pr)
    
    kk <- seq_along(map$all.bet)
    names(kk) <- as.character(map$all.bet)
    k <- 0
    for(i in 1:nmodel){
      id.a <- alp.index.lo(map, i)
      if(is.null(id.a)){
        next
      }
      
      id.b <- map$bet[[i]]
      for(j in id.a){
        k <- k + 1
        tmp <- g.alp.lam[, k]
        for(l in id.b){
          h[j, l] <- -g.alp.bet[[foo(j, l)]] + 
            sum(tmp * g.bet.lam[, kk[as.character(l)]] * pr2)
          h[l, j] <- t(h[j, l])
        }
      }
    }
    
    rm(g.alp.bet)
    gc()
    
    k <- 0
    for(j in map$all.alp){
      k <- k + 1
      tmp <- g.alp.lam[, k]
      for(l in map$all.bet){
        
        if(!is.na(h[j, l])){
          next
        }
        
        h[j, l] <- sum(tmp * g.bet.lam[, kk[as.character(l)]] * pr2)
        h[l, j] <- t(h[j, l])
      }
    }
  }
  
  ####### ell.bet.bet
  
  g.bet2 <- gfunction.bet2.lo(para, map, ref, delta, lam, pr)
  kk <- seq_along(map$all.bet)
  names(kk) <- as.character(map$all.bet)
  k <- 0
  for(i in 1:nmodel){
    id.b <- map$bet[[i]]
    
    for(j in id.b){
      k <- k + 1
      tmp <- g.bet.lam[, k]
      for(l in id.b){
        if(l < j){
          next
        }
        
        h[j, l] <- -g.bet2[[foo(j, l)]] + 
          sum(tmp * g.bet.lam[, kk[as.character(l)]] * pr2)
        h[l, j] <- t(h[j, l])
      }
    }
  }
  
  rm(g.bet2)
  gc()
  
  k <- 0
  for(j in map$all.bet){
    k <- k + 1
    tmp <- g.bet.lam[, k]
    for(l in map$all.bet){
      if(l <= j){
        next
      }
      
      if(!is.na(h[j, l])){
        next
      }
      
      h[j, l] <- sum(tmp * g.bet.lam[, kk[as.character(l)]] * pr2)
      h[l, j] <- t(h[j, l])
    }
  }
  
  h[map$all.bet, map$all.bet] <- h[map$all.bet, map$all.bet] - inv.V
  
  ##########
  
  colnames(h) <- names(para)
  rownames(h) <- names(para)
  
  h
  
}
