

reorganize.lm <- function(fit, map){
  
  score <- fit$score
  para <- fit$coefficients
  vcov <- fit$vcov
  
  sigma2 <- para[map$the[1]]
  para <- para[map$the[-1]]
  vcov <- vcov[map$the[-1], map$the[-1]]
  rownames(vcov) <- names(para)
  colnames(vcov) <- names(para)
  
  fit <- list(coefficients = para, vcov = vcov, sigma2 = sigma2, score = score)
  fit
  
}
