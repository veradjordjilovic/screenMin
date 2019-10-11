require(mvtnorm)
require(Rsolnp)
source("CoreScreenMin.R")

# September 21, 2019


################################################
#  Example
#
#  Set up the parameters 
d <- 3.25   # signal-to-noise ratio
p0 <- 0.85  # proportion of (0,0) hypotheses (H_i such that both H_i1 and H_i2 are true)
p1 <- 0.1   # proportion of (1,0) hypotheses (H_i such that exactly one of H_i1 and H_i2 is true)
p2 <- 0.05  # proportion of (1,1) hypotheses (H_i such that both H_i1 and H_i2 are false)
level <- 0.05  # level at which FWER is to be controlled
mh <- 100   # the number of hypotheses

# Find c* that maximizes the probability of rejecting a false hypothesis
# under the costraint that the FWER is bounded by alpha
threshold  <- solnp(c(0.01), f,  ineqfun = ineq, mu = d, pi0 = p0, pi1 = p1, pi2 = p2, alpha = level, m = mh,
              ineqLB = c(-Inf), ineqUB = c(0), LB=10^(-10), UB =0.05)

threshold$pars
threshold$values


# Comparison with the threshold found by exhausting the FWER
threshold_FWER <- uniroot(safe_t_exact, c(0.001, 0.07), mu = d,  pi0 = p0, pi1 = p1, pi2 = p2, 
        alpha = level,
        m = mh, tol = 10^(-9)) 

threshold_FWER$root 
# they are the same

#############################################
# Plot power and FWER bound as a function of the threshold

f.V <- Vectorize(f, vectorize.args = "x")

# Function to compute FWER estimate (similar to safe_t_exact)
f.t <- function (x, mu, pi0, pi1, pi2, alpha, m){
  # FWER of ScreenMin
  # c = x
  
  # F(c)
  t <- Prob(c=x, mu = mu) 
  
  # E(S)
  expS <- (pi2*(2*t-t^2)+ pi0*(2*x-x^2) + pi1*(t+x-x*t))*m 
  
  # F(alpha/E(S))
  u <- Prob(c = alpha/expS, mu = mu) 
  
  if ((alpha/expS<=x)) { Prej <-alpha/expS*u/ (x + t - x*t) } else{Prej <- (x*u+alpha/expS*t -x*t)/(x + t - x*t)}
  res <- 1-(1-Prej)^expS
  return(res)
}
f.t.V <- Vectorize(f.t, vectorize.args = "x")


x <- c(seq(from=0.0001, to= threshold$pars, length.out = 100),
       seq(from= threshold$pars, to= 0.05, length.out = 100))

y <- -f.V(x, mu = d, 
          pi0 = p0, pi1 = p1, pi2 = p2, alpha = level, m = mh)

z <- f.t.V(x, mu = d, pi0 = p0, pi1 = p1, pi2 = p2, alpha = level, m = mh)

# Power of the Bonferroni procedure for comparison
Bonf.power<- (Prob(c = level/mh, mu = d))^2

plot(x,z, xlim= range(x), ylim =  c(min(c(y,z, Bonf.power)),1.1* max(c(y,z))), type="l",
     xlab = "c",
     ylab= "", lty = 2, bty="l", lwd = 2)
points(x,y, type = "l", lty = 1, lwd = 2)
abline(h = level, lty=3, col = "gray", lwd=2)
abline(v = threshold_FWER$root, lty = 3, col = "gray", lwd=2)
abline(h = Bonf.power, lty = 6, lwd = 2)
legend("topright", lty=c(1, 6, 2), lwd = 2, bty="n", 
       legend= c("ScreenMin power",   "Bonferroni power", "ScreenMin FWER"))




