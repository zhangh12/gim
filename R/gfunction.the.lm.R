

gfunction.the.lm <- function(para, map, data){
  
  data$'(Intercept)' <- 1
  
  nmodel <- length(map$bet)
  
  sigma <- para[map$the[1]]
  the <- para[map$the[-1]]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  
  g.the <- list()
  
  nthe <- length(the)
  n <- nrow(data)
  
  nlam <- max(map$lam)
  
  gt <- matrix(0, nrow = n, ncol = nlam)
  offset <- max(map$the)
  for(i in 1:nmodel){
    id.tau <- map$alp[[i]][1]
    gt[, id.tau - offset] <- 1
  }
  
  g.the[[1]] <- gt
  rm(gt)
  
  for(j in 1:nthe){
    gt <- matrix(0, nrow = n, ncol = nlam)
    fx0 <- fx[, names(the)[j]]
    for(i in 1:nmodel){
      id.tau <- map$alp[[i]][1]
      tau <- para[id.tau]
      
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
      
      rx <- as.matrix(data[, names(gam), drop = FALSE])
      
      delta <- as.vector(fx %*% the - rx %*% gam)
      
      gt[, id.tau - offset] <- 2 * delta * fx0
      id <- c(id.a, id.b)
      gt[, id - offset] <- rx[, names(para)[id], drop = FALSE] * fx0
      
      rm(id.tau, id.a, id.b, alp.exist)
    }
    g.the[[j+1]] <- gt
    rm(gt)
  }
  
  g.the
  
}


