

gfunction.bet.lm <- function(para, map, ref){
  
  nmodel <- length(map$bet)
  
  sigma <- para[map$the[1]]
  the <- para[map$the[-1]]
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  
  g.bet <- list()
  
  n <- nrow(ref)
  nlam <- max(map$lam)
  offset <- max(map$the)
  
  k <- 0
  for(i in 1:nmodel){
    id.a <- alp.index.lm(map, i)
    alp.exist <- !is.null(id.a)
    if(alp.exist){
      alp <- para[id.a]
    }else{
      alp <- NULL
    }
    
    id.b <- map$bet[[i]]
    bet <- para[id.b]
    gam <- c(alp, bet)
    
    rx <- as.matrix(ref[, names(gam), drop = FALSE])
    
    delta <- as.vector(fx %*% the - rx %*% gam)
    
    for(j in id.b){
      rx0 <- rx[, names(para)[j]]
      gb <- matrix(0, nrow = n, ncol = nlam)
      if(alp.exist){
        gb[, id.a - offset] <- -rx[, names(alp), drop = FALSE] * rx0
      }
      gb[, id.b - offset] <- -rx[, names(bet), drop = FALSE] * rx0
      k <- k + 1
      g.bet[[k]] <- gb
      rm(gb)
    }
  }
  
  g.bet
  
}


