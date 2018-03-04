

gfunction.the.lm <- function(para, para.id, data){
  
  data$'(Intercept)' <- 1
  
  id.lam <- para.id$id.lam
  id.the <- para.id$id.the
  id.alp <- para.id$id.alp
  id.bet <- para.id$id.bet
  
  nmodel <- nrow(id.bet)
  
  sigma <- para[id.the$start[1]]
  the <- para[(id.the$start[1]+1):id.the$end[1]]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  
  g.the <- list()
  
  nthe <- length(the)
  n <- nrow(data)
  
  nlam <- max(id.lam)
  
  gt <- matrix(0, nrow = n, ncol = nlam)
  offset <- max(id.the)
  for(i in 1:nmodel){
    id.tau <- id.alp$start[i]
    gt[, id.tau - offset] <- 1
  }
  
  g.the[[1]] <- gt
  rm(gt)
  
  for(j in 1:nthe){
    gt <- matrix(0, nrow = n, ncol = nlam)
    fx0 <- fx[, names(the)[j]]
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

