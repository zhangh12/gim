
# g.sigma = 0, as the first component
gfunction.the.lm <- function(para, map, ref){
  
  nmodel <- length(map$bet)
  
  the <- para[map$the[-1]]
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  
  nthe <- length(the)
  n <- nrow(ref)
  
  nlam <- max(map$lam)
  
  offset <- max(map$the)
  
  g.the <- list()
  g.the[[1]] <- matrix(0, nrow = n, ncol = nlam)
  
  for(j in 1:nthe){
    gt <- matrix(0, nrow = n, ncol = nlam)
    
    fx0 <- fx[, names(the)[j]]
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
      
      id <- c(id.a, id.b)
      gt[, id - offset] <- rx[, names(para)[id], drop = FALSE] * fx0
      
    }
    g.the[[j + 1]] <- gt
    rm(gt)
  }
  
  g.the
  
}


