## September 26, 2019
# Simulation study


require(mvtnorm)
require(Rsolnp)
require(MultiMed)
source("AdaFilterBonf.R")


Prob <- function(c=0.05, mu = 0){
  # Probability that a (one-sided) p-value of a normal test statistic is below c, when the true
  # signal to noise ratio is mu
  pnorm(qnorm(1-c), mean= mu, lower.tail=F)
}

f <- function (x, mu, pi0, pi1, pi2, alpha, m){
  # Probability of rejecting a false hypothesis by ScreenMin
  #
  # Input:
  #   x... value between 0 and 1 to be used as screening threshold 
  #   mu... positive value representing the SNR of the false hypotheses
  #   pi0... proportion of H_i such that both H_i1 and H_i2 are true
  #   pi1... proportion of H_i such that exactly one of H_i1 and H_i2 is true
  #   pi2... proportion of H_i such that both H_i1 and H_i2 are false
  #   alpha... level at which FWER is to be controlled
  #   m...  number of hypotheses
  
  # 
  # Output:
  #   Minus the estimated probability of rejecting a false hypothesis when the selection threshold is x.
  #
  t <- Prob(c = x, mu = mu)
  expS <- (pi2*(2*t-t^2)+ pi0*(2*x-x^2) + pi1*(t+x-x*t))*m
  u <- Prob(c = alpha/expS, mu = mu)
  if (alpha < x*expS) {
    res <- u^2 } else{ res <- t^2 + 2*t*(u-t)
    }
  return(-res) }



safe_t_exact <- function (x, mu, pi0, pi1, pi2, alpha, m){
  # Evaluate FWER for a given selection threshold
  # Input:
  #   See f 
  #
  # Output:
  #  A logical value stating whether the FWER estimate is below alpha or not
  t <- Prob(c = x, mu = mu)
  expS <- (pi2*(2*t-t^2)+ pi0*(2*x-x^2) + pi1*(t+x-x*t))*m
  u <- Prob(c = alpha/expS, mu = mu)
  if (alpha/expS<=x) {res <- alpha*u/(expS*(t+x-t*x))} else{
    res <- (alpha*t/expS + x*u-x*t)/(t+x-x*t) 
  }
  res <- 1 - (1-res)^expS - alpha
  return(res)}

ineq <- function (x, mu, pi0, pi1, pi2, alpha, m){
  y <- safe_t_exact(x, mu, pi0, pi1, pi2, alpha, m)
  return(c(y))}


SimFWER <- function (mu, pi0, pi1, pi2, alpha, m, B = 100, fac=1, sseed = 123, rho = 0){
  #  Simulation study to assess FWER and power for ScreenMin
  #
  #  Args:
  #   mu... positive value representing the SNR of the false hypotheses
  #   pi0... proportion of H_i such that both H_i1 and H_i2 are true
  #   pi1... proportion of H_i such that exactly one of H_i1 and H_i2 is true
  #   pi2... proportion of H_i such that both H_i1 and H_i2 are false
  #   alpha... level at which FWER is to be controlled
  #   m...  number of hypotheses
  #   B... the number of Monte Carlo runs
  #   fac... the ratio between the SNR in the two columns of the p-values
  #
  #  Output:
  #    A vector of size 4 giving computed FWER, FDR and associated powers. 
  #
  
  # prepare thresholds
  
  threshold  <- solnp(c(0.01), f,  ineqfun = ineq, mu = mu, pi0 = pi0, pi1 = pi1,
                      pi2 = pi2, alpha = alpha, m = m,
                      ineqLB = c(-Inf), ineqUB = c(0), LB=10^(-10), UB =0.05)
  
  opt.threshold <- threshold$pars
  def.threshold <- alpha/m
  
  
  
  # true parameters
  m1 <- round(m*(pi1+pi2))
  m2 <- round(m*pi2)
  mu1 <- c(rep(mu, m1), rep(0, m-m1))
  mu2 <- c(rep(mu*fac, m2), rep(0, m-m2)) # unequal power
  sigma <- matrix(rho, ncol=m, nrow=m) + diag(rep(1-rho, m))
  
  
  # true and false hypotheses
  TH <- mu1*mu2 == 0
  FH <- !TH  
  
  
  
  # results matrix
  FWER<- matrix(0, ncol = 5, nrow = B)
  power <- matrix(0, ncol = 5, nrow = B)
  colnames(FWER) <- colnames(power) <- c("Optimal threshold SM", "Default threshold SM",
                                         "Bonf.",
                                         "MCP_S", "AdaFilter")
  
  set.seed(sseed)
  for (i in 1:B){
    if (i%%100 ==0) print(i)
    if (rho == 0) {
      p1 <-  pnorm(rnorm(m, mean = mu1, sd= rep(1,m)), lower.tail=F)
      p2 <- pnorm(rnorm(m, mean = mu2, sd= rep(1,m)), lower.tail=F)} else {
        p1 <- pnorm(rmvnorm(1, mean = mu1, sigma= sigma), lower.tail=F)
        p2 <- pnorm(rmvnorm(1, mean = mu2, sigma= sigma), lower.tail=F)}
    
    minp <- pmin(p1, p2)
    maxp <- pmax(p1, p2)
    
    # screenMin optimal
    S <- which(minp <=  opt.threshold)
    SM <- which(maxp < alpha/ length(S))
    SM <- intersect(S, SM)
    psmo <- rep(1, m)
    names(psmo) <- 1:m
    psmo[SM] <- maxp[SM]*length(SM)
    
    # screenMin default
    Sd <- which(minp <= def.threshold)
    SMd<- which(maxp <= alpha/ length(Sd))
    SMd <- intersect(Sd, SMd)
    psmd <- rep(1, m)
    names(psmd) <- 1:m
    psmd[SMd] <- maxp[SMd]*length(SMd)
    
    # Bonferroni
    pb <- as.numeric(pmin(maxp*m, 1))
    names(pb) <- 1:m
    
    # Heller 2018
    ph <- medTest.SBMH(p1, p2, MCP.type="FWER", t1=alpha/2, t2=alpha/2, lambda=0)
    
    # AdaFilter
    
    Ada.p <- AdaFilterBonf(p1, p2, alpha=alpha)
    
    # results
    pval <- cbind(psmo, psmd, pb, ph, Ada.p)
    rej <- apply( pval, 2, function(x) x <= alpha)
    FWER[i, ] <- apply(rej, 2, function(x) sum(x[TH])>0)
    if (sum(FH) ==0) {power[i, ] <- rep(0, 5)} else{
      power[i, ] <- apply(rej, 2, function(x) sum(x[FH])/ sum(FH))}
  }
  
  FWER <- colMeans(FWER)
  power <- colMeans(power)
  
  return(list("FWER" = FWER, "Power" = power))
} 


Rcpp::sourceCpp(code='
          #include <RcppArmadillo.h>
          // [[Rcpp::depends(RcppArmadillo)]]
          
          using namespace Rcpp;
          
          // [[Rcpp::export]]
          arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
          int ncols = sigma.n_cols;
          arma::mat Y = arma::randn(n, ncols);
          return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
          }'
)

SimFWER_hd <- function (mu, pi0, pi1, pi2, alpha, m, B = 100, fac=1, sseed = 123, rho = 0){
  #  Computes FWER and power for ScreenMin when m is large
  #
  #  Args:
  #   mu... positive value representing the SNR of the false hypotheses
  #   pi0... proportion of H_i such that both H_i1 and H_i2 are true
  #   pi1... proportion of H_i such that exactly one of H_i1 and H_i2 is true
  #   pi2... proportion of H_i such that both H_i1 and H_i2 are false
  #   alpha... level at which FWER is to be controlled
  #   m...  number of hypotheses
  #   B... the number of Monte Carlo runs
  #   fac... the ratio between the SNR in the two columns of the p-values
  #
  #  Output:
  #    A vector of size 4 giving computed FWER, FDR and associated powers. 
  #
  
  # prepare thresholds
  
  threshold  <- solnp(c(0.01), f,  ineqfun = ineq, mu = mu, pi0 = pi0, pi1 = pi1,
                      pi2 = pi2, alpha = alpha, m = m,
                      ineqLB = c(-Inf), ineqUB = c(0), LB=10^(-10), UB =0.05)
  
  opt.threshold <- threshold$pars
  def.threshold <- alpha/m
  
  
  
  
  m1 <- round(m*(pi1+pi2))
  m2 <- round(m*pi2)
  mu1 <- c(rep(mu, m1), rep(0, m-m1))
  mu2 <- c(rep(mu*fac, m2), rep(0, m-m2)) # unequal power
  sigma <- matrix(rho, ncol=m, nrow=m) + diag(rep(1-rho, m))
  
  # true and false hypotheses
  TH <- mu1*mu2 == 0
  FH <- !TH
  
  
  # results matrix
  FWER<- matrix(0, ncol = 5, nrow = B)
  power <- matrix(0, ncol = 5, nrow = B)
  colnames(FWER) <- colnames(power) <- c("Optimal threshold SM", "Default threshold SM",
                                         "Bonf.",
                                         "MCP_S", "AdaFilter")
  set.seed(sseed)
  P1 = mvrnormArma(B, mu1, sigma)
  P2 = mvrnormArma(B, mu2, sigma)
  P1 <- pnorm(P1, lower.tail=F)
  P2 <- pnorm(P2, lower.tail=F)
  
  for (i in 1:B){
    if (i%%100 ==0) print(i)
    
    p1 <- P1[i, ]
    p2 <- P2[i, ]
    minp <- pmin(p1, p2)
    maxp <- pmax(p1, p2)
    # screenMin optimal
    S <- which(minp <=  opt.threshold)
    SM <- which(maxp < alpha/ length(S))
    SM <- intersect(S, SM)
    psmo <- rep(1, m)
    names(psmo) <- 1:m
    psmo[SM] <- maxp[SM]*length(SM)
    
    # screenMin default
    Sd <- which(minp <= def.threshold)
    SMd<- which(maxp <= alpha/ length(Sd))
    SMd <- intersect(Sd, SMd)
    psmd <- rep(1, m)
    names(psmd) <- 1:m
    psmd[SMd] <- maxp[SMd]*length(SMd)
    
    # Bonferroni
    pb <- as.numeric(pmin(maxp*m, 1))
    names(pb) <- 1:m
    
    # Heller 2018
    ph <- medTest.SBMH(p1, p2, MCP.type="FWER", t1=alpha/2, t2=alpha/2, lambda=0)
    
    # AdaFilter
    
    Ada.p <- AdaFilterBonf(p1, p2, alpha=alpha)
    
    # results
    pval <- cbind(psmo, psmd, pb, ph, Ada.p)
    rej <- apply( pval, 2, function(x) x <= alpha)
    FWER[i, ] <- apply(rej, 2, function(x) sum(x[TH])>0)
    if (sum(FH) ==0) {power[i, ] <- rep(0, 5)} else{
      power[i, ] <- apply(rej, 2, function(x) sum(x[FH])/ sum(FH))}
    
  }
  
  FWER <- colMeans(FWER)
  power <- colMeans(power)
  
  return(list("FWER" = FWER, "Power" = power))
}  



#res1 <- SimFWER(mu = 3, pi0 = 0.85, pi1 = 0.1, pi2 = 0.05, alpha = 0.05,
#               m = 50, B = 10000, fac=1)
#res1
