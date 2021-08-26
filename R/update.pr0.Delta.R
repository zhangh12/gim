

update.pr0.Delta <- function(para, map, ref, Delta, ncase, nctrl, pr0){
  
  if(is.null(pr0) + is.null(Delta) == 1){
    stop('Message for Kai (LoL): pr0 and Delta should be NULL (or not NULL) at the same time')
  }
  
  if(is.null(pr0) + is.null(Delta) == 0){
    return(list(pr0 = pr0, Delta = Delta))
  }
  
  tilt <- tilt.cc(para, map, ref)
  Delta <- tilt$Delta
  delta <- tilt$delta
  
  lam <- para[map$lam]
  g <- gfunction.cc(para, map, ref, Delta, delta, ncase, nctrl)
  n <- nrow(g)
  pr0 <- as.vector(1 / (1 + g %*% lam)) / n
  
  return(list(pr0 = pr0, Delta = Delta))
  
}

