

gfunction.bet.lm <- function(para, para.id, data){
  
  data$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.bet)
  
  sigma <- para[id.the$start[1]]
  the <- para[(id.the$start[1]+1):id.the$end[1]]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  
  g.bet <- list()
  
  n <- nrow(data)
  nlam <- max(id.lam)
  offset <- max(id.the)
  
  k <- 0
  for(i in 1:nmodel){
    id.tau <- id.alp$start[i]
    tau <- para[id.tau]
    
    id.a <- alp.index.lm(id.alp, i)
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
    
    delta <- as.vector(fx %*% the - rx %*% gam)
    
    for(j in id.b){
      rx0 <- rx[, names(para)[j]]
      gb <- matrix(0, nrow = n, ncol = nlam)
      gb[, id.tau - offset] <- -2 * delta * rx0
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

