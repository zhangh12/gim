
check.neg.def <- function(h){
  
  ei <- eigen(h)
  v <- ei$values
  if(all(v < 0)){
    return(h)
  }
  
  message('Levenberg–Marquardt')
  v[v >= 0] <- max(v[v < 0])
  
  U <- ei$vectors
  
  U %*% diag(v) %*% t(U)
  
}
