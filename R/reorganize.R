

reorganize <- function(fit, para.id){
  
  para <- fit$coefficients
  vcov <- fit$vcov
  
  id.the <- para.id$id.the
  sigma2 <- para[min(id.the)]
  para <- para[(min(id.the)+1):max(id.the)]
  vcov <- vcov[(min(id.the)+1):max(id.the), (min(id.the)+1):max(id.the)]
  rownames(vcov) <- names(para)
  colnames(vcov) <- names(para)
  
  fit <- list(coefficients = para, vcov = vcov, sigma2 = sigma2)
  fit
  
}
