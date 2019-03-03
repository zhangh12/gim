

gfunction.the2.lm <- function(para, map, ref, lam){
  
  nmodel <- length(map$bet)
  
  g.the2 <- list()
  
  the <- para[map$the]
  nthe <- length(map$the)
  n <- nrow(ref)
  
  nlam <- max(map$lam)
  offset <- max(map$the)
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  for(j in 2:nthe){
    fxj <- ref[, names(the)[j]]
    for(l in j:nthe){
      fxl <- ref[, names(the)[l]]
      gt <- matrix(0, nrow = n, ncol = nlam)
      gt[, map$all.tau - offset] <- 2 * fxj * fxl
      
      g.the2[[foo(j,l)]] <- gt %*% lam
      rm(gt)
    }
  }
  
  g.the2
  
}


