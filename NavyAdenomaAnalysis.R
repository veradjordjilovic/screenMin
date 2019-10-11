library(MultiMed)
source("AdaFilterBonf.R")

data(NavyAdenoma)
colnames(NavyAdenoma)[1:5]
colnames(NavyAdenoma)[c(6:9,154)]
colnames(NavyAdenoma)[155]
plot(density(NavyAdenoma$Fish))
hist(NavyAdenoma$Fish)

# prevalence of adenoma in a population of interest
prev <- 0.228
p <- sum(NavyAdenoma$Adenoma==1)/nrow(NavyAdenoma)
p

# weights for cases and controls
w <- rep(NA, nrow(NavyAdenoma))
w[NavyAdenoma$Adenoma == 1] <- prev/p
w[NavyAdenoma$Adenoma == 0] <- (1-prev)/(1-p)
table(w)


# Exploratory analysis 
boxplot(BMI ~ Adenoma, data = NavyAdenoma,names= c("Controls", "Cases"), main="BMI", frame.plot=T, 
        border=c( "black","red"))
boxplot(Fish ~ Adenoma, data= NavyAdenoma ,names= c("Controls", "Cases"), main="Fish intake", frame.plot=T, 
        border=c( "black","red"))
boxplot(Age ~ Adenoma, data= NavyAdenoma ,names= c("Controls", "Cases"), main="Age", frame.plot=T, 
        border=c( "black","red"))


# Total effect of fish consumption on the risk of adenoma is weak
fit_exposure <- glm( Adenoma ~ Fish + BMI + Female + Age + Smoking, data= NavyAdenoma)
summary(fit_exposure)

# P-values for the association with fish intake (p1) and the risk of adenoma (p2) for each metabolite
p1 <- p2 <- numeric()
for (i in 6:154){
  fitM <- lm( NavyAdenoma[, i] ~ as.matrix(NavyAdenoma[, 1:5]), weights = w )
  p1[i] <- summary(fitM)$coefficients[2, 4]  
  res <- residuals(fitM)
  fitY <- glm(Adenoma ~ as.matrix(NavyAdenoma[, c(1:5, i)]), 
              family="binomial", data=NavyAdenoma)
  p2[i] <- summary(fitY)$coefficients[7, 4]
}

names(p1) <- 1:length(p1)
names(p2) <- 1:length(p2)
p1 <- p1[-c(1:5)]
p2 <- p2[-c(1:5)]

# ScreenMin procedure
minp <- pmin(p1, p2)
maxp <- pmax(p1, p2)
S <- names(which(pmin(p1,p2) < 0.05/149))
length(S)
p.adjust(maxp[S], method="bonferroni")
SM <- names(which(pmax(p1, p2) < 0.05/ length(S)))
SM <- intersect(S, SM)

# AdaFilter procedure
adaf <- AdaFilterBonf(p1, p2, alpha=0.05)
min(adaf)
which(adaf == min(adaf))
names(p1)[which(adaf == min(adaf))]

# Adaptive threshold of AdaFilter
alpha <- 0.05
sortmin <- sort(minp)
temp <- sortmin *(1:length(p1))
ind.temp <- which(temp>alpha)
if (length(ind.temp)>0) {k.prim <- ind.temp[1]
if (sortmin[k.prim]*(k.prim-1) <= alpha) {k <- k.prim}else{
  k <- k.prim-1}} else {
    k <- m}
adapt.threshold <- alpha/k
adapt.S <- names(which(pmin(p1,p2) <= adapt.threshold))
length(adapt.S)

# MultiMed
alpha <- 0.05
multimedtest <- medTest.SBMH(p1, p2, MCP.type="FWER", t1=alpha/2, t2=alpha/2, lambda=0)
min(multimedtest)
which(multimedtest == min(multimedtest))
names(p1)[which(multimedtest == min(multimedtest))]


# Would any threshold have lead to some rejections?
SM <- numeric()
sortmin <- sort(minp)
for (i in (1:length(p1))){
  c <- sortmin[i]
  S <- names(which(pmin(p1,p2) <= c))
  temp <- names(which(pmax(p1, p2) <= 0.05/ length(S)))
  SM[i] <- length(intersect(S, temp))
}

# Which thresholds lead to rejections?
sortmin[which(SM >0)]

# Compare with default screenmin
alpha/149

# Why don't we find any mediators?
require(globaltest)
gtest <- gt(NavyAdenoma$Adenoma, NavyAdenoma[, 6:149], NavyAdenoma[, 1:5], model="logistic")
# The global test suggests that the group of 149 metabolites is asscoiated with the risk of adenoma

# Plot the scatterplot of p-values
plot("n", bty="l", main= "",
     xlab=expression(p[1]), ylab=expression(p[2]), xlim=c(0,1), ylim=c(0,1))

shadedX = c(-0.01,alpha, alpha, -0.01)
shadedY = c(-0.01,-0.01,1,1)
polygon(shadedX,shadedY,col="bisque2", border = NA) 
shadedY = c(-0.01,alpha, alpha, -0.01)
shadedX = c(0,0,1,1)
polygon(shadedX,shadedY,col="bisque2", border = NA) 
shadedX = c(-0.01,alpha, alpha, -0.01)
shadedY = c(-0.01, -0.01, alpha, alpha)
polygon(shadedX,shadedY,col="bisque4", border = NA) 
points(p1, p2, pch=20)

# There seems to be little evidence that the fish intake has impact on the metabolites
