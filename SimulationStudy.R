# September 26 2019
# Simulation study complete

source("SimulationStudyFunctions.R")

######################################################################
# Independence + Figure 1 
######################################################################

# setting 1
r1 <- SimFWER(mu = 3, pi0 = 0.95, pi1 = 0, pi2 = 0.05, alpha = 0.05,
              m = 200, B = 1000, fac=1)
r2 <- SimFWER(mu = 3, pi0 = 0.85, pi1 = 0.1, pi2 = 0.05, alpha = 0.05,
              m = 200, B = 1000, fac=1)
r3 <- SimFWER(mu = 3, pi0 = 0.75, pi1 = 0.2, pi2 = 0.05, alpha = 0.05,
              m = 200, B = 1000, fac=1)
r4 <- SimFWER(mu = 3, pi0 = 0.65, pi1 = 0.3, pi2 = 0.05, alpha = 0.05,
              m = 200, B = 1000, fac=1)
r5 <- SimFWER(mu = 3, pi0 = 0.55, pi1 = 0.4, pi2 = 0.05, alpha = 0.05,
              m = 200, B = 1000, fac=1)


# setting 2: m = 10000
mh <- 10000

r6 <- SimFWER(mu = 4, pi0 = 0.95, pi1 = 0, pi2 = 0.05, alpha = 0.05,
              m = mh, B = 1000, fac=1)
r7 <- SimFWER(mu = 4, pi0 = 0.85, pi1 = 0.1, pi2 = 0.05, alpha = 0.05,
              m = mh, B = 1000, fac=1)
r8 <- SimFWER(mu = 4, pi0 = 0.75, pi1 = 0.2, pi2 = 0.05, alpha = 0.05,
              m = mh, B = 1000, fac=1)
r9 <- SimFWER(mu = 4, pi0 = 0.65, pi1 = 0.3, pi2 = 0.05, alpha = 0.05,
              m = mh, B = 1000, fac=1)
r10 <- SimFWER(mu = 4, pi0 = 0.55, pi1 = 0.4, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=1)

# setting 3: m=200 and unequal SNR for the two studies
mh <- 200
fc <- 2


r11 <- SimFWER(mu = 3, pi0 = 0.95, pi1 = 0, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc)
r12 <- SimFWER(mu = 3, pi0 = 0.85, pi1 = 0.1, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc)
r13 <- SimFWER(mu = 3, pi0 = 0.75, pi1 = 0.2, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc)
r14 <- SimFWER(mu = 3, pi0 = 0.65, pi1 = 0.3, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc)
r15 <- SimFWER(mu = 3, pi0 = 0.55, pi1 = 0.4, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc)


Power1 <- rbind(r1$Power, r2$Power, r3$Power, r4$Power,r5$Power)
Power2 <- rbind(r6$Power, r7$Power, r8$Power, r9$Power,r10$Power)
Power3 <- rbind(r11$Power, r12$Power, r13$Power, r14$Power,r15$Power)



FWER1 <- rbind(r1$FWER, r2$FWER, r3$FWER, r4$FWER, r5$FWER)
FWER2 <- rbind(r6$FWER, r7$FWER, r8$FWER, r9$FWER, r10$FWER)
FWER3 <- rbind(r11$FWER, r12$FWER, r13$FWER, r14$FWER,r15$FWER)


xax <- c(0,0.1, 0.2, 0.3, 0.4)
t <-1
#pdf ("Figures/SSIndepUpdated.pdf", height=6.6, width=10)
par(mfrow=c(2,3))
plot(y= FWER1[,1], x=xax, xlab=expression(pi[1]),
     ylab="FWER",ylim=c(0,0.05), bty="l",
     type="b", pch=22, lty=2, cex=t,
     main= expression(paste(italic(m)==200, ", equal SNR")))
points(y= FWER1[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= FWER1[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= FWER1[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= FWER1[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

legend("topleft", legend=c("Oracle SM", "AdaFilter",
                           "Default SM", "MCP_S",  "Bonf."),
       pch=c(22,3, 17,21, 23), bty="n")

plot(y= FWER3[,1], x=xax, xlab=expression(pi[1]),
     ylab="FWER",ylim=c(0,0.05), bty="l",
     type="b", pch=22, lty=2, cex=t,
     main= expression(paste(italic(m)==200, ", unequal SNR")))
points(y= FWER3[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= FWER3[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= FWER3[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= FWER3[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

plot(y= FWER2[,1], x=xax, xlab=expression(pi[1]),
     ylab="FWER",ylim=c(0,0.05), bty="l",
     type="b", pch=22, lty=2, cex=t,
     main= expression(paste(italic(m)==10000, ", equal SNR")))
points(y= FWER2[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= FWER2[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= FWER2[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= FWER2[,5], x=xax, type="b",  lty=2, pch=3, cex=t)


plot(y= Power1[,1], x=xax, xlab=expression(pi[1]),
     ylab="Power",ylim=c(0,0.5), bty="l",
     type="b", pch=22, lty=2, cex=t)
points(y= Power1[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= Power1[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= Power1[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= Power1[,5], x=xax, type="b",  lty=2, pch=3, cex=t)


plot(y= Power3[,1], x=xax, xlab=expression(pi[1]),
     ylab="Power",ylim=c(0,0.7), bty="l",
     type="b", pch=22, lty=2, cex=t)
points(y= Power3[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= Power3[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= Power3[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= Power3[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

plot(y= Power2[,1], x=xax, xlab=expression(pi[1]),
     ylab="Power",ylim=c(0,0.45), bty="l",
     type="b", pch=22, lty=2, cex=t)
points(y= Power2[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= Power2[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= Power2[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= Power2[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

#dev.off()
#########################################################
# Dependence + Figure 2
#########################################################

# setting 4:  m=200, rho=0.3
r <- 0.3
r16 <- SimFWER(mu = 3, pi0 = 0.95, pi1 = 0, pi2 = 0.05, alpha = 0.05,
               m = 200, B = 1000, fac=1, rho = r)
r17 <- SimFWER(mu = 3, pi0 = 0.85, pi1 = 0.1, pi2 = 0.05, alpha = 0.05,
               m = 200, B = 1000, fac=1,  rho = r)
r18 <- SimFWER(mu = 3, pi0 = 0.75, pi1 = 0.2, pi2 = 0.05, alpha = 0.05,
               m = 200, B = 1000, fac=1, rho = r)
r19 <- SimFWER(mu = 3, pi0 = 0.65, pi1 = 0.3, pi2 = 0.05, alpha = 0.05,
               m = 200, B = 1000, fac=1, rho =r)
r20 <- SimFWER(mu = 3, pi0 = 0.55, pi1 = 0.4, pi2 = 0.05, alpha = 0.05,
               m = 200, B = 1000, fac=1, rho = r)

# setting 5: m = 200 and rho = 0.8
mh <- 200
fc <- 1
r <- 0.8
r21 <- SimFWER(mu = 3, pi0 = 0.95, pi1 = 0, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc, rho = r)
r22 <- SimFWER(mu = 3, pi0 = 0.85, pi1 = 0.1, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc, rho = r)
r23 <- SimFWER(mu = 3, pi0 = 0.75, pi1 = 0.2, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc, rho = r)
r24 <- SimFWER(mu = 3, pi0 = 0.65, pi1 = 0.3, pi2= 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc, rho = r)
r25 <- SimFWER(mu = 3, pi0 = 0.55, pi1 = 0.4, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc, rho = r)

# setting 6: m = 200 and rho = 0.8, fc = 2 
mh <- 200
fc <- 2
r <- 0.8 

r26 <- SimFWER(mu = 3, pi0 = 0.95, pi1 = 0, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc, rho = r)
r27 <- SimFWER(mu = 3, pi0 = 0.85, pi1 = 0.1, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc, rho = r)
r28 <- SimFWER(mu = 3, pi0 = 0.75, pi1 = 0.2, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc, rho = r)
r29 <- SimFWER(mu = 3, pi0 = 0.65, pi1 = 0.3, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc, rho = r)
r30 <- SimFWER(mu = 3, pi0 = 0.55, pi1 = 0.4, pi2 = 0.05, alpha = 0.05,
               m = mh, B = 1000, fac=fc, rho = r)


Power4 <- rbind(r16$Power, r17$Power, r18$Power, r19$Power, r20$Power)
Power5 <- rbind(r21$Power, r22$Power, r23$Power, r24$Power, r25$Power)
Power6 <- rbind(r26$Power, r27$Power, r28$Power, r29$Power, r30$Power)


FWER4 <- rbind(r16$FWER, r17$FWER, r18$FWER, r19$FWER, r20$FWER)
FWER5 <- rbind(r21$FWER, r22$FWER, r23$FWER, r24$FWER, r25$FWER)
FWER6 <- rbind(r26$FWER, r27$FWER, r28$FWER, r29$FWER, r30$FWER)


xax <- c(0,0.1, 0.2, 0.3, 0.4)
t <-1
#pdf ("Figures/SSdepUpdated.pdf", height=6.6, width=10)
par(mfrow=c(2,3))
plot(y= FWER4[,1], x=xax, xlab=expression(pi[1]),
     ylab="FWER",ylim=c(0,0.05), bty="l",
     type="b", pch=22, lty=2, cex=t,
     main= expression(paste(rho==0.3, ", equal SNR")))
points(y= FWER4[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= FWER4[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= FWER4[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= FWER4[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

legend("topleft", legend=c("Oracle SM", "AdaFilter",
                           "Default SM", "MCP_S",  "Bonf."),
       pch=c(22,3, 17,21, 23), bty="n")

plot(y= FWER5[,1], x=xax, xlab=expression(pi[1]),
     ylab="FWER",ylim=c(0,0.05), bty="l",
     type="b", pch=22, lty=2, cex=t,
     main= expression(paste(rho==0.8, ", equal SNR")))
points(y= FWER5[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= FWER5[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= FWER5[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= FWER5[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

plot(y= FWER6[,1], x=xax, xlab=expression(pi[1]),
     ylab="FWER",ylim=c(0,0.05), bty="l",
     type="b", pch=22, lty=2, cex=t,
     main= expression(paste(rho==0.8, ", unequal SNR")))
points(y= FWER6[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= FWER6[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= FWER6[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= FWER6[,5], x=xax, type="b",  lty=2, pch=3, cex=t)


plot(y= Power4[,1], x=xax, xlab=expression(pi[1]),
     ylab="Power",ylim=1.1*c(0,max(Power4)), bty="l",
     type="b", pch=22, lty=2, cex=t)
points(y= Power4[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= Power4[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= Power4[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= Power4[,5], x=xax, type="b",  lty=2, pch=3, cex=t)


plot(y= Power5[,1], x=xax, xlab=expression(pi[1]),
     ylab="Power",ylim=1.1*c(0,max(Power5)), bty="l",
     type="b", pch=22, lty=2, cex=t)
points(y= Power5[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= Power5[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= Power5[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= Power5[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

plot(y= Power6[,1], x=xax, xlab=expression(pi[1]),
     ylab="Power",ylim = 1.1*c(0,max(Power6)), bty="l",
     type="b", pch=22, lty=2, cex=t)
points(y= Power6[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= Power6[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= Power6[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= Power6[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

#dev.off()
######################################################################
# Dependence high-dimensional 
######################################################################

# setting 7: high dimensional setting m=10000, rho = 0.3 (can take hours to run)
mh <- 10000
fc <- 1
r <- 0.3 
d <- 4

r31 <- SimFWER_hd(mu = d, pi0 = 0.95, pi1 = 0, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 10, fac=fc, rho = r)
r32 <- SimFWER_hd(mu = d, pi0 = 0.85, pi1 = 0.1, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)
r33 <- SimFWER_hd(mu = d, pi0 = 0.75, pi1 = 0.2, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)
r34 <- SimFWER_hd(mu = d, pi0 = 0.65, pi1 = 0.3, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)
r35 <- SimFWER_hd(mu = 3, pi0 = 0.55, pi1 = 0.4, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)

# setting 8: high dimensional setting m=10000, rho = 0.8
mh <- 10000
fc <- 1
r <- 0.8
d <- 4

r36 <- SimFWER_hd(mu = d, pi0 = 0.95, pi1 = 0, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)
r37 <- SimFWER_hd(mu = d, pi0 = 0.85, pi1 = 0.1, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)
r38 <- SimFWER_hd(mu = d, pi0 = 0.75, pi1 = 0.2, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)
r39 <- SimFWER_hd(mu = d, pi0 = 0.65, pi1 = 0.3, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)
r40 <- SimFWER_hd(mu = 3, pi0 = 0.55, pi1 = 0.4, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)

# setting 9: high dimensional setting m=10000, rho = 0.8
mh <- 10000
fc <- 2
r <- 0.8
d <- 4

r41 <- SimFWER_hd(mu = d, pi0 = 0.95, pi1 = 0, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)
r42 <- SimFWER_hd(mu = d, pi0 = 0.85, pi1 = 0.1, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)
r43 <- SimFWER_hd(mu = d, pi0 = 0.75, pi1 = 0.2, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)
r44 <- SimFWER_hd(mu = d, pi0 = 0.65, pi1 = 0.3, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)
r45 <- SimFWER_hd(mu = 3, pi0 = 0.55, pi1 = 0.4, pi2 = 0.05, alpha = 0.05,
                  m = mh, B = 1000, fac=fc, rho = r)

# Figure 
Power7 <- rbind(r31$Power, r32$Power, r33$Power, r34$Power, r35$Power)
Power8 <- rbind(r36$Power, r37$Power, r38$Power, r39$Power, r40$Power)
Power9 <- rbind(r41$Power, r42$Power, r43$Power, r44$Power, r45$Power)


FWER7 <- rbind(r31$FWER, r32$FWER, r33$FWER, r34$FWER, r35$FWER)
FWER8 <- rbind(r36$FWER, r37$FWER, r38$FWER, r39$FWER, r40$FWER)
FWER9 <- rbind(r41$FWER, r42$FWER, r43$FWER, r44$FWER, r45$FWER)


xax <- c(0,0.1, 0.2, 0.3, 0.4)
t <-1
#pdf ("Figures/SShdUpdated.pdf", height=6.6, width=10)
par(mfrow=c(2,3))
plot(y= FWER7[,1], x=xax, xlab=expression(pi[1]),
     ylab="FWER",ylim=c(0,0.05), bty="l",
     type="b", pch=22, lty=2, cex=t,
     main= expression(paste(rho==0.3, ", equal SNR")))
points(y= FWER7[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= FWER7[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= FWER7[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= FWER7[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

legend("topleft", legend=c("Oracle SM", "AdaFilter",
                           "Default SM", "MCP_S",  "Bonf."),
       pch=c(22,3, 17,21, 23), bty="n")

plot(y= FWER8[,1], x=xax, xlab=expression(pi[1]),
     ylab="FWER",ylim=c(0,0.05), bty="l",
     type="b", pch=22, lty=2, cex=t,
     main= expression(paste(rho==0.8, ", equal SNR")))
points(y= FWER8[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= FWER8[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= FWER8[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= FWER8[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

plot(y= FWER9[,1], x=xax, xlab=expression(pi[1]),
     ylab="FWER",ylim=c(0,0.05), bty="l",
     type="b", pch=22, lty=2, cex=t,
     main= expression(paste(rho==0.8, ", unequal SNR")))
points(y= FWER9[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= FWER9[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= FWER9[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= FWER9[,5], x=xax, type="b",  lty=2, pch=3, cex=t)


plot(y= Power7[,1], x=xax, xlab=expression(pi[1]),
     ylab="Power",ylim=1.1*c(0,max(Power7)), bty="l",
     type="b", pch=22, lty=2, cex=t)
points(y= Power7[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= Power7[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= Power7[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= Power7[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

#legend("topright", legend=c("Oracle SM", "AdaFilter",
#                           "Default SM", "MCP_S",  "Bonf."),
#       pch=c(22,3, 17,21, 23), bty="n")

plot(y= Power8[,1], x=xax, xlab=expression(pi[1]),
     ylab="Power",ylim=1.1*c(0,max(Power8)), bty="l",
     type="b", pch=22, lty=2, cex=t)
points(y= Power8[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= Power8[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= Power8[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= Power8[,5], x=xax, type="b",  lty=2, pch=3, cex=t)
hhy
plot(y= Power9[,1], x=xax, xlab=expression(pi[1]),
     ylab="Power",ylim = 1.1*c(0,max(Power9)), bty="l",
     type="b", pch=22, lty=2, cex=t)
points(y= Power9[,2], x=xax, type="b", lty=2, pch=17, cex=t)
points(y= Power9[,3], x=xax, type="b", lty=2, pch=23, cex=t, bg="grey")
points(y= Power9[,4], x=xax, type="b",  lty=2, pch=21, cex=t)
points(y= Power9[,5], x=xax, type="b",  lty=2, pch=3, cex=t)

#dev.off()

#save.image(file="SSUpdated26092019.RData")
