# 1  October 2019

ScreenMinBonf <- function(p1, p2, c){
  # SreenMin procedure
  #
  # Args:
  #  p1... an m-length vector of p-values (study1)
  #  p2... second m-length vector of p-values (study2)
  #  c... threshold for selection
  #
  # Returns:
  # an m-length vector of Bonferroni adjusted p-values
  #
  minp <- pmin(p1, p2)
  maxp <- pmax(p1, p2)
  
  m <- length(p1)
   
  S <- which(minp <=  c)
  SM <- which(maxp <= alpha/ length(S))
  SM <- intersect(S, SM)
  adj.pval <- rep(1, m)
  names(adj.pval) <- 1:m
  adj.pval[SM] <- maxp[SM]*length(SM) 
  return(adj.pval)
}

# Example:
m <- 100
alpha <- 0.05
p1 <-  pnorm(rnorm(m, mean = 2, sd= rep(1,m)), lower.tail=F)
p2 <-  pnorm(rnorm(m, mean = 3, sd= rep(1,m)), lower.tail=F)

p.val.def <- ScreenMinBonf(p1, p2, alpha/m)
which(p.val.def <= alpha)
