

gfunction.the.alp.cc <- function(para, map, ref, Delta, delta, ncase, nctrl, xi, pr){
  
  nmodel <- length(map$bet)
  
  the <- para[map$the]
  
  g.the.alp <- list()
  
  nthe <- length(the)
  n <- nrow(ref)
  
  const <- list()
  for(i in 1:nmodel){
    id <- c(alp.index.cc(map, i), map$bet[[i]])
    gam <- para[id]
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    
    rho.i <- ncase[i, i] / nctrl[i, i]
    const[[i]] <- -rx * (Delta * rho.i * delta[, i] / (1 + rho.i * delta[, i])^2)
  }
  
  nlam <- max(map$lam)
  offset <- max(map$the)
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  for(j in 1:nthe){
    fxj <- ref[, names(the)[j]]
    for(i in 1:nmodel){
      id.a <- alp.index.cc(map, i)
      if(is.null(id.a)){
        next
      }
      id <- c(alp.index.cc(map, i), map$bet[[i]])
      tmp <- const[[i]] * fxj
      for(l in id.a){
        fxl <- ref[, names(para)[l]]
        gt <- matrix(0, nrow = n, ncol = nlam - 1)
        gt[, id - offset] <- tmp * fxl
        #g.the.alp[[foo(j,l)]] <- gt
        g.the.alp[[foo(j,l)]] <- t(gt %*% xi) %*% pr
        rm(gt)
      }
    }
  }
  
  if(length(g.the.alp) == 0){
    g.the.alp <- NULL
  }
  
  g.the.alp
  
}


