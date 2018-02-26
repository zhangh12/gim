

reorganize.lo <- function(fit, para.id){
  
  para <- fit$coefficients
  vcov <- fit$vcov
  
  id.the <- para.id$id.the
  para <- para[min(id.the):max(id.the)]
  vcov <- vcov[min(id.the):max(id.the), min(id.the):max(id.the)]
  rownames(vcov) <- names(para)
  colnames(vcov) <- names(para)
  
  fit <- list(coefficients = para, vcov = vcov)
  fit
  
}
