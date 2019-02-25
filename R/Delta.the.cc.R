

Delta.the.cc <- function(para, map, ref, Delta){
  
  the <- para[map$the]
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  
  fx * Delta
  
}
