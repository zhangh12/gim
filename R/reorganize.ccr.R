

reorganize.ccr <- function(fit, map){
  
  score <- fit$score
  para <- fit$coefficients
  vcov <- fit$vcov
  
  id <- c(map$ome, map$the)
  para <- para[id]
  vcov <- vcov[id, id]
  rownames(vcov) <- names(para)
  colnames(vcov) <- names(para)
  
  np <- length(para)
  A <- cbind(0, diag(1, np - 1))
  A[1, 1] <- 1
  
  para1 <- as.vector(A %*% para)
  names(para1) <- names(para)[-1]
  vcov1 <- A %*% vcov %*% t(A)
  rownames(vcov1) <- names(para1)
  colnames(vcov1) <- names(para1)
  
  fit <- list(coefficients = para1, vcov = vcov1, score = score)
  fit
  
}
