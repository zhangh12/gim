
# data is not necessary if numerical derivative is not used
hess.cc <- function(para, map, data, ref, inv.V, bet0, sample.info, outcome){
  
  if(0){
    h <- numDeriv::jacobian(score.cc, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, sample.info = sample.info, outcome = outcome)
    colnames(h) <- names(para)
    rownames(h) <- names(para)
    return(h)
  }
  
  np <- length(para)
  h <- matrix(NA, np, np)
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  n <- nrow(data)
  
  lam <- para[map$lam]
  xi <- lam[-1]
  
  tilt <- tilt.cc(para, map, ref)
  Delta <- tilt$Delta
  delta <- tilt$delta
  
  ncase <- sample.info$ncase
  nctrl <- sample.info$nctrl
  
  g <- gfunction.cc(para, map, ref, Delta, delta, ncase, nctrl)
  g1 <- g[, -1, drop = FALSE]
  g.the <- gfunction.the.cc(para, map, ref, Delta, delta, ncase, nctrl)
  g.alp <- gfunction.alp.cc(para, map, ref, Delta, delta, ncase, nctrl)
  g.bet <- gfunction.bet.cc(para, map, ref, Delta, delta, ncase, nctrl)
  
  g.xi <- as.vector(g1 %*% xi)
  g.the.xi <- gfunction.the.xi.cc(g.the, xi)
  g.alp.xi <- gfunction.alp.xi.cc(g.alp, xi)
  g.bet.xi <- gfunction.alp.xi.cc(g.bet, xi)
  
  pr <- as.vector(1/(1+g %*% lam))
  pr2 <- pr^2
  
  ####### ell.lam.lam
  
  h[map$lam, map$lam] <- t(g) %*% (g * pr2)
  
  ####### ell.lam[1].the
  
  Delta.the <- Delta.the.cc(para, map, ref, Delta)
  nthe <- length(map$the)
  
  h[map$the, map$lam[1]] <- -t(Delta.the) %*% ((1 + g.xi) * pr2) + t(g.the.xi) %*% (g[, 1] * pr2)
  h[map$lam[1], map$the] <- t(h[map$the, map$lam[1]])
  
  ####### ell.xi.the
  
  h[map$lam[-1], map$the] <- t(g1) %*% ((lam[1] * Delta.the + g.the.xi) * pr2)
  for(i in seq_along(map$the)){
    h[map$lam[-1], map$the[i]] <- h[map$lam[-1], map$the[i]] - t(g.the[[i]]) %*% pr
  }
  h[map$the, map$lam[-1]] <- t(h[map$lam[-1], map$the])
  
  ####### ell.lam[1].alp
  
  if(!is.null(map$all.alp)){
    h[map$all.alp, map$lam[1]] <- t(g.alp.xi) %*% (g[, 1] * pr2)
    h[map$lam[1], map$all.alp] <- t(h[map$all.alp, map$lam[1]])
  }
  
  ####### ell.lam[1].bet
  
  h[map$all.bet, map$lam[1]] <- t(g.bet.xi) %*% (g[, 1] * pr2)
  h[map$lam[1], map$all.bet] <- t(h[map$all.bet, map$lam[1]])
  
  ####### ell.xi.alp
  
  if(!is.null(map$all.alp)){
    h[map$lam[-1], map$all.alp] <- t(g1) %*% (g.alp.xi * pr2)
    nalp <- length(map$all.alp)
    for(i in seq_along(map$all.alp)){
      h[map$lam[-1], map$all.alp[i]] <- h[map$lam[-1], map$all.alp[i]] - t(g.alp[[i]]) %*% pr
    }
    
    h[map$all.alp, map$lam[-1]] <- t(h[map$lam[-1], map$all.alp])
  }
  
  ####### ell.xi.bet
  
  h[map$lam[-1], map$all.bet] <- t(g1) %*% (g.bet.xi * pr2)
  nbet <- length(map$all.bet)
  for(i in seq_along(map$all.bet)){
    h[map$lam[-1], map$all.bet[i]] <- h[map$lam[-1], map$all.bet[i]] - t(g.bet[[i]]) %*% pr
  }
  
  h[map$all.bet, map$lam[-1]] <- t(h[map$lam[-1], map$all.bet])
  
  ####### ell.the.the
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  g.the2 <- gfunction.the2.cc(para, map, ref, Delta, delta, ncase, nctrl, xi, pr)
  the <- para[map$the]
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  for(j in seq_along(map$the)){
    id.j <- map$the[j]
    fxj <- fx[, names(the)[j]]
    tmp1 <- Delta * fxj
    tmp2 <- tmp1 * pr
    for(l in j:nthe){
      id.l <- map$the[l]
      fxl <- fx[, names(the)[l]]
      h[id.j, id.l] <- -lam[1] * sum(tmp2 * fxl) -
        g.the2[[foo(j,l)]] +
        sum((lam[1] * tmp1 + as.vector(g.the.xi[, j])) * 
              (lam[1] * Delta * fxl + as.vector(g.the.xi[, l])) * pr2)
      h[id.l, id.j] <- t(h[id.j, id.l])
    }
  }
  rm(g.the2)
  gc()
  
  ####### ell.the.alp
  
  if(!is.null(map$all.alp)){
    g.the.alp <- gfunction.the.alp.cc(para, map, ref, Delta, delta, ncase, nctrl, xi, pr)
    for(j in seq_along(map$the)){
      id.j <- map$the[j]
      fxj <- fx[, names(the)[j]]
      tmp <- as.vector(lam[1] * Delta * fxj + g.the.xi[, j]) * pr2
      k <- 0
      for(i in 1:nmodel){
        id.a <- alp.index.cc(map, i)
        if(is.null(id.a)){
          next
        }
        
        for(l in id.a){
          k <- k + 1
          h[id.j, l] <- -g.the.alp[[foo(j, l)]] + 
            sum(g.alp.xi[, k] * tmp)
          h[l, id.j] <- t(h[id.j, l])
        }
      }
    }
    
    rm(g.the.alp)
    gc()
  }
  
  ####### ell.the.bet
  
  g.the.bet <- gfunction.the.bet.cc(para, map, ref, Delta, delta, ncase, nctrl, xi, pr)
  for(j in seq_along(map$the)){
    id.j <- map$the[j]
    fxj <- fx[, names(the)[j]]
    tmp <- as.vector(lam[1] * Delta * fxj + g.the.xi[, j]) * pr2
    k <- 0
    for(i in 1:nmodel){
      id.b <- map$bet[[i]]
      
      for(l in id.b){
        k <- k + 1
        h[id.j, l] <- -g.the.bet[[foo(j, l)]] + 
          sum(g.bet.xi[, k] * tmp)
        h[l, id.j] <- t(h[id.j, l])
      }
    }
  }
  
  rm(g.the.bet)
  gc()
  
  ####### ell.alp.alp
  
  if(!is.null(map$all.alp)){
    g.alp2 <- gfunction.alp2.cc(para, map, ref, Delta, delta, ncase, nctrl, xi, pr)
    
    kk <- seq_along(map$all.alp)
    names(kk) <- as.character(map$all.alp)
    k <- 0
    for(i in 1:nmodel){
      id.a <- alp.index.cc(map, i)
      if(is.null(id.a)){
        next
      }
      
      for(j in id.a){
        k <- k + 1
        tmp <- g.alp.xi[, k]
        for(l in id.a){
          if(l < j){
            next
          }
          
          h[j, l] <- -g.alp2[[foo(j, l)]] + 
            sum(tmp * g.alp.xi[, kk[as.character(l)]] * pr2)
          h[l, j] <- t(h[j, l])
        }
      }
    }
    
    rm(g.alp2)
    gc()
    
    k <- 0
    for(j in map$all.alp){
      k <- k + 1
      tmp <- g.alp.xi[, k]
      for(l in map$all.alp){
        if(l <= j){
          next
        }
        
        if(!is.na(h[j, l])){
          next
        }
        
        h[j, l] <- sum(tmp * g.alp.xi[, kk[as.character(l)]] * pr2)
        h[l, j] <- t(h[j, l])
      }
    }
  }
  
  ####### ell.alp.bet
  
  if(!is.null(map$all.alp)){
    g.alp.bet <- gfunction.alp.bet.cc(para, map, ref, Delta, delta, ncase, nctrl, xi, pr)
    
    kk <- seq_along(map$all.bet)
    names(kk) <- as.character(map$all.bet)
    k <- 0
    for(i in 1:nmodel){
      id.a <- alp.index.cc(map, i)
      if(is.null(id.a)){
        next
      }
      
      id.b <- map$bet[[i]]
      for(j in id.a){
        k <- k + 1
        tmp <- g.alp.xi[, k]
        for(l in id.b){
          h[j, l] <- -g.alp.bet[[foo(j, l)]] + 
            sum(tmp * g.bet.xi[, kk[as.character(l)]] * pr2)
          h[l, j] <- t(h[j, l])
        }
      }
    }
    
    rm(g.alp.bet)
    gc()
    
    k <- 0
    for(j in map$all.alp){
      k <- k + 1
      tmp <- g.alp.xi[, k]
      for(l in map$all.bet){
        
        if(!is.na(h[j, l])){
          next
        }
        
        h[j, l] <- sum(tmp * g.bet.xi[, kk[as.character(l)]] * pr2)
        h[l, j] <- t(h[j, l])
      }
    }
  }
  
  ####### ell.bet.bet
  
  g.bet2 <- gfunction.bet2.cc(para, map, ref, Delta, delta, ncase, nctrl, xi, pr)
  kk <- seq_along(map$all.bet)
  names(kk) <- as.character(map$all.bet)
  k <- 0
  for(i in 1:nmodel){
    id.b <- map$bet[[i]]
    
    for(j in id.b){
      k <- k + 1
      tmp <- g.bet.xi[, k]
      for(l in id.b){
        if(l < j){
          next
        }
        
        h[j, l] <- -g.bet2[[foo(j, l)]] + 
          sum(tmp * g.bet.xi[, kk[as.character(l)]] * pr2)
        h[l, j] <- t(h[j, l])
      }
    }
  }
  
  rm(g.bet2)
  gc()
  
  k <- 0
  for(j in map$all.bet){
    k <- k + 1
    tmp <- g.bet.xi[, k]
    for(l in map$all.bet){
      if(l <= j){
        next
      }
      
      if(!is.na(h[j, l])){
        next
      }
      
      h[j, l] <- sum(tmp * g.bet.xi[, kk[as.character(l)]] * pr2)
      h[l, j] <- t(h[j, l])
    }
  }
  
  h[map$all.bet, map$all.bet] <- h[map$all.bet, map$all.bet] - inv.V
  
  ####################
  
  colnames(h) <- names(para)
  rownames(h) <- names(para)
  h
  
}

