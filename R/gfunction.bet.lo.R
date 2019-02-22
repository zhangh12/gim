

gfunction.bet.lo <- function(para, map, data){
  
  data$'(Intercept)' <- 1
  
  nmodel <- length(map$bet)
  
  the <- para[map$the]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  
  g.bet <- list()
  
  n <- nrow(data)
  nlam <- max(map$lam)
  offset <- max(map$the)
  
  k <- 0
  for(i in 1:nmodel){
    
    id.a <- alp.index.lo(map, i)
    alp.exist <- !is.null(id.a)
    if(alp.exist){
      alp <- para[id.a]
    }else{
      alp <- NULL
    }
    
    id.b <- map$bet[[i]]
    bet <- para[id.b]
    gam <- c(alp, bet)
    
    rx <- as.matrix(data[, names(gam), drop = FALSE])
    tmp <- as.vector(exp(rx %*% gam))
    yi <- tmp/(1+tmp)
    
    for(j in id.b){
      rx0 <- rx[, names(para)[j]]
      gb <- matrix(0, nrow = n, ncol = nlam)
      
      if(alp.exist){
        gb[, id.a - offset] <- -rx[, names(alp), drop = FALSE] * (rx0 * yi * (1-yi))
      }
      gb[, id.b - offset] <- -rx[, names(bet), drop = FALSE] * (rx0 * yi * (1-yi))
      k <- k + 1
      g.bet[[k]] <- gb
      rm(gb)
    }
  }
  
  g.bet
  
}


