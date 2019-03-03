

gfunction.tau.lm <- function(para, map, ref){
  
  nmodel <- length(map$bet)
  n <- nrow(ref)
  
  g.tau <- list()
  
  nlam <- max(map$lam)
  offset <- max(map$the)
  
  for(i in 1:nmodel){
    id.tau <- map$alp[[i]][1]
    tau <- para[id.tau]
    
    gt <- matrix(0, nrow = n, ncol = nlam)
    gt[, id.tau - offset] <- -1
    g.tau[[i]] <- gt
    rm(gt)
  }
  
  g.tau
  
}


