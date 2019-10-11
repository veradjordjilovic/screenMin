# Core functions for the ScreenMin threshold analysis

# October 11, 2019

Prob <- function(c=0.05, mu = 0){
  # Probability that a (one-sided) p-value of a normal test statistic is below c, when the true
  # signal to noise ratio is mu
  pnorm(qnorm(1-c), mean= mu, lower.tail=F)
}

f <- function (x, mu, pi0, pi1, pi2, alpha, m){
  # Probability of rejecting a false hypothesis by ScreenMin
  #
  # Input:
  #   x... value between 0 and 1 representing the screening threshold 
  #   mu... positive value representing the SNR of the false hypotheses
  #   pi0... proportion of H_i such that both H_i1 and H_i2 are true
  #   pi1... proportion of H_i such that exactly one of H_i1 and H_i2 is true
  #   pi2... proportion of H_i such that both H_i1 and H_i2 are false
  #   alpha... level at which FWER is to be controlled
  #   m...  number of hypotheses
  
  # 
  # Output:
  #   Minus the approximated probability of rejecting a false hypothesis when the threshold is x.
  #
  t <- Prob(c = x, mu = mu)
  expS <- (pi2*(2*t-t^2)+ pi0*(2*x-x^2) + pi1*(t+x-x*t))*m
  u <- Prob(c = alpha/expS, mu = mu)
  if (alpha < x*expS) {
    res <- u^2 } else{ res <- t^2 + 2*t*(u-t)
    }
  return(-res) }

safe_t_exact <- function (x, mu, pi0, pi1, pi2, alpha, m){
  # Evaluate whether the unconditional FWER of ScreenMin is below alpha for a given threshold c
  # Input:
  #   See f 
  #
  # Output:
  #  Difference between the (upper bound of) ScreenMin FWER and alpha
  #
  # F(c)
  t <- Prob(c = x, mu = mu)
  
  # E(S)
  expS <- (pi2*(2*t-t^2)+ pi0*(2*x-x^2) + pi1*(t+x-x*t))*m
  
  # F(alpha/E(S))
  u <- Prob(c = alpha/expS, mu = mu)
  
  # P_0(alpha/E(S), c)
  if (alpha/expS<=x) {res <- alpha*u/(expS*(t+x-t*x))} else{
    res <- (alpha*t/expS + x*u-x*t)/(t+x-x*t) 
  }
  
  # ScreenMin FWER - alpha
  res <- 1 - (1-res)^expS - alpha
  return(res)}

ineq <- function (x, mu, pi0, pi1, pi2, alpha, m){
  # Auxuliary function for optimization
  y <- safe_t_exact(x, mu, pi0, pi1, pi2, alpha, m)
  return(c(y))}
