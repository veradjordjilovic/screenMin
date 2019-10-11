# 22 September 2019

AdaFilterBonf <- function(p1, p2, alpha){
  # AdaFilter procedure of Wang et al.(2018)
  #
  # Args:
  #  p1... an m-length vector of p-values (study1)
  #  p2... second m-length vector of p-values (study2)
  #  alpha... level at which FWER is to be controlled
  #
  # Returns:
  # an m-length vector of Bonferroni adjusted p-values
  #
  minp <- pmin(p1, p2)
  maxp <- pmax(p1, p2)
  
  m <- length(p1)
  
  sortmin <- sort(minp)
  temp <- sortmin *(1:m)
  ind.temp <- which(temp>alpha)
  if (length(ind.temp)>0) {k.prim <- ind.temp[1]
                           if (sortmin[k.prim]*(k.prim-1) <= alpha) {k <- k.prim}else{
                             k <- k.prim-1}} else {
                             k <- m}
  threshold <- alpha/k
  pb <- as.numeric(pmin(maxp*k, 1))
  names(pb) <- 1:m
  
  return(pb)
}

# Example:
m <- 100
alpha <- 0.05
p1 <-  pnorm(rnorm(m, mean = 2, sd= rep(1,m)), lower.tail=F)
p2 <-  pnorm(rnorm(m, mean = 3, sd= rep(1,m)), lower.tail=F)
p.val.ada <- AdaFilterBonf(p1, p2, alpha)
which(p.val.ada<=alpha)
