

reorganize.cml <- function(fit, map){
  
  score <- fit$score
  para <- fit$coefficients
  vcov <- fit$vcov
  
  para <- para[map$the]
  vcov <- vcov[map$the, map$the]
  rownames(vcov) <- names(para)
  colnames(vcov) <- names(para)
  
  fit <- list(coefficients = para, vcov = vcov, score = score)
  fit
  
}
