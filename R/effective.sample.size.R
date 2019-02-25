
# effective sample size for case-control data
# harmonic mean of (n0, n1)
effective.sample.size <- function(n0, n1){
  
  2/(1/n0 + 1/n1)
  
}
