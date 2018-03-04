

gfunction.alp.lo <- function(para, para.id, data){
  
  data$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.bet)
  
  the <- para[id.the$start[1]:id.the$end[1]]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  
  g.alp <- list()
  
  n <- nrow(data)
  nlam <- max(id.lam)
  offset <- max(id.the)
  
  k <- 0
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
    
    rx <- as.matrix(data[, names(gam), drop = FALSE])
    tmp <- as.vector(exp(rx %*% gam))
    yi <- tmp/(1+tmp)
    
    if(!alp.exist){
      rm(id.a, id.b, alp.exist)
      next
    }
    
    for(j in id.a){
      rx0 <- rx[, names(para)[j]]
      ga <- matrix(0, nrow = n, ncol = nlam)
      
      ga[, id.a - offset] <- -rx[, names(alp), drop = FALSE] * (rx0 * yi * (1-yi))
      ga[, id.b - offset] <- -rx[, names(bet), drop = FALSE] * (rx0 * yi * (1-yi))
      k <- k + 1
      g.alp[[k]] <- ga
      rm(ga)
    }
    
    rm(id.a, id.b, alp.exist)
  }
  
  if(length(g.alp) == 0){
    g.alp <- NULL
  }
  
  g.alp
  
}

