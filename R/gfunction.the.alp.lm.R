
# sigma is not in 'the'
gfunction.the.alp.lm <- function(para, map, ref, lam){
  
  nmodel <- length(map$bet)
  
  g.the.alp <- list()
  
  the <- para[map$the]
  nthe <- length(the)
  offset <- max(map$the)
  n <- nrow(ref)
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  zero <- rep(0, n)
  for(j in 2:nthe){
    fxj <- ref[, names(the)[j]]
    for(i in 1:nmodel){
      id.tau <- map$alp[[i]][1]
      #g.the.alp[[foo(j, id.tau)]] <- zero
      
      id.a <- alp.index.lm(map, i)
      if(is.null(id.a)){
        next
      }
      
      for(l in id.a){
        fxl <- ref[, names(para)[l]]
        g.the.alp[[foo(j,l)]] <- -2 * fxj * fxl * lam[id.tau - offset]
      }
    }
  }
  
  g.the.alp
  
}


