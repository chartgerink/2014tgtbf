setwd("D:/Dropbox/projects/2014tgtbf/figures")

# Figure 1
N <- 10000
i <- seq(1, N, 1)
fcv <- qf(p=1-.05 ,df1=1,df2=97)

# Null
beta <- pf(q=fcv, df1=1, df2=97)
quadpoints <- (i/(N+1))*beta
fquad <- qf(quadpoints, df1=1, df2=97)
temp <- 1-pf(fquad, df1=1, df2=97)
pquad0 <- (temp-.05)/(1-.05)

# Small
beta <- pf(q=fcv, df1=1, df2=97, ncp=(.01/(1-.01))*(100+1+1))
N <- 10000
i <- seq(1, N, 1)
quadpoints <- (i/(N+1))*beta
fquad <- qf(quadpoints, df1=1, df2=97, ncp=(.01/(1-.01))*(100+1+1))
temp <- 1-pf(fquad, df1=1, df2=97)
pquadS <- (temp-.05)/(1-.05)

# Medium
beta <- pf(q=fcv, df1=1, df2=97, ncp=(.0625/(1-.0625))*(100+1+1))
quadpoints <- (i/(N+1))*beta
fquad <- qf(quadpoints, df1=1, df2=97, ncp=(.0625/(1-.0625))*(100+1+1))
temp <- 1-pf(fquad, df1=1, df2=97)
pquadM <- (temp-.05)/(1-.05)

# Large
beta <- pf(q=fcv, df1=1, df2=97, ncp=(.14/(1-.14))*(100+1+1))
quadpoints <- (i/(N+1))*beta
fquad <- qf(quadpoints, df1=1, df2=97, ncp=(.14/(1-.14))*(100+1+1))
temp <- 1-pf(fquad, df1=1, df2=97)
pquadL <- (temp-.05)/(1-.05)

pdf('Fig1.pdf', width=7, height=8)
par(mai=c(1,1,.2,.2))
plot(density(pquadL, kernel="gaussian", bw="SJ", adjust=1), xlim=c(0.01,1), xaxs='i', lwd=4,
     frame.plot=T, 
     main="",
     #      ylim=c(0,4.5),
     xaxs="i",
     yaxs="i",
     xlab="P-value",
     ylab = "Density",
     cex.axis=.8,
     cex.lab=1,
     col = "black", las=1, bty = 'L')
lines(density(pquadM, kernel="gaussian", bw="SJ", adjust=1), xlim=c(0.01,1), xaxs='i', ylim=c(0,10), lty=3, lwd=3)
lines(density(pquadS, kernel="gaussian", bw="SJ",adjust=.5), xlim=c(0.01,1), xaxs='i', ylim=c(0,10), lwd=2)
abline(h=1)
legend(x=.6,y=5,legend=c(latex2exp('Large effect, $\\eta^2=.14$'),
                         latex2exp('Medium effect, $\\eta^2=.06$'),
                         latex2exp('Small effect, $\\eta^2=.01$'),
                         latex2exp('No effect, $\\eta^2=0$')),
       cex=1, lty=c(1,3,1,1),
       col = "black", lwd=c(4,3,2,1), box.lwd=0 , bty='n')
dev.off()