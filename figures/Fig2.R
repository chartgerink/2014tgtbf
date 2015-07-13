# Adapted from code provided by Michèle Nuijten

setwd("D:/Dropbox/projects/2014tgtbf/figures")
library(ggplot2)
library(latex2exp)

# pdf("Fig2.pdf", width = 7, height = 4)
setEPS()
postscript("whatever.eps") # change file name
par(mar = c(1, 2, 3, 4))
plot(dnorm(seq(-2, 2, .001), 0, .5), type = 'l',
     bty = 'n', xaxt = 'n', yaxt = 'n')
lines(dnorm(seq(-2, 2, .001), 1, .5))

dev.off



par(mar = c(0, 0, 0, 0))
par(xpd = FALSE)

# distribution H0
h0 <- rnorm(n = 1e7, 0, 1)
density.h0 <- density(h0)
q40.h0 <- quantile(h0, .4)
q95.h0 <- quantile(h0, .95)
q100.h0 <- quantile(h0, 1)

dd.h0 <- with(density.h0, data.frame(x, y))

# distribution h1
h1 <- rnorm(n = 1e7, 2, 1)
density.h1 <- density(h1)
q0.h1 <- quantile(h1, 0)
q35.h1 <- quantile(h1, 0.375)

dd.h1 <- with(density.h1, data.frame(x, y))

# Plot h0
ggplot(dd.h0, aes(x, y)) +
  geom_line() + 
  geom_line(data = dd.h1, aes(x, y)) + 
  geom_ribbon(data = subset(dd.h0, x > q95.h0 & x < q100.h0), 
              aes(ymax = y), ymin = 0,
              fill = "black", colour = NA, alpha = 0.5) +
  geom_ribbon(data = subset(dd.h1, x > q0.h1 & x < q35.h1), 
              aes(ymax = y), ymin = 0,
              fill = "red", colour = NA, alpha = 0.5) + 
  geom_ribbon(data = subset(dd.h0, x > q40.h0 & x < q100.h0), 
              aes(ymax = y), ymin = 0,
              fill = "grey", colour = NA, alpha = 0.5) + 
  theme_bw()
  annotate("text", label = latex2exp('\\beta'),
            x = 0, y = .1)
annotate("text", label = "plot mpg vs. wt", x = 2, y = 15, size = 8, colour = "red")



  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank(),
    axis.line = element_blank()
  )


# plot distribution under H0
plot(density(H0),
     xlim = c(-4, 6 + d), 
     ylim = c(0.015, .5),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "", main = "")

# distribution H1
H1 <- rt(n=1e7,df=df,ncp=ncp)

# plot distribution under H0
lines(density(H1))

# calculate critical value
cv <- qt(.95,df = df)

# shade power area
draws <- rt(2e7,df = df,ncp = ncp)
dens <- density(draws)

max <- max(H0 + ncp)

x1 <- min(which(dens$x <= cv))  
x2 <- max(which(dens$x <  max))

with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), 
                   col=adjustcolor("gray",alpha.f=0.5), border=NA))

text(2.6,.2,expression(1-beta))
#----------------------------------

# shade alpha area
draws <- rt(2e7,df=df,ncp=0)
dens <- density(draws)

max <- max(H0)

x1 <- min(which(dens$x >= cv))  
x2 <- max(which(dens$x <  max))

with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), 
                   col=adjustcolor("darkgrey",alpha.f=0.5), border=NA))

text(cv+.25,.03,expression(alpha))

#----------------------------------

# shade publication bias area
densy <- pub*density(H1)$y
densx <- density(H1)$x
min <- min(densx)

x1 <- min(which(densx >= min))  
x2 <- max(which(densx <  cv))

polygon(x=c(densx[c(x1,x1:x2,x2)]), y= c(0, densy[x1:x2], 0), 
        col=adjustcolor("grey20",alpha.f=0.5), border=NA)


lines(densx[x1:x2],pub*density(H1)$y[x1:x2],lty=4)

#----------------------------------

# determine estimated effect sizes

# range
points <- 1000
p <- seq(1,points)
p <- p/(points+1)

t <- qt(p,df,ncp)
tcv <- cv

# select t values conditional on rejection of H0
# calculate weight
# calculate expected t
t_sig <- t[t>tcv]
weight_t_sig <- 1/points
exp_t_sig <- sum(t_sig*weight_t_sig)
D1 <- exp_t_sig/(weight_t_sig*length(t_sig))

# select t values conditional on acceptation of H0
# calculate weight
# calculate expected t
t_ns <- t[t<=tcv]
weight_t_ns <- pub*(1/points)
exp_t_ns <- sum(t_ns*weight_t_ns)
D0 <- exp_t_ns/(weight_t_ns*length(t_ns))

# expected t value in literature
# will get overestimated when less n.s. findings are published
t_lit <- (exp_t_sig+exp_t_ns)/(weight_t_sig*length(t_sig)+weight_t_ns*length(t_ns))
D <- t_lit

#----------------------------------

# (re)draw lines

# plot distribution under Ha
lines(density(H1))

# redo line distribution H0
lines(density(H0))

#-----------

abline(v=cv)
mtext("cv",3,at=cv,line=.5)

#-----------

# indicate mean of H0
arrows(mean(H0),max(dt(H0,df=df)),mean(H0),-.1,length=0,lty=2)
mtext("0",1,at=mean(H0),line=.5)

# indicate mean of H1
arrows(ncp,max(dt(H1,df=df)),ncp,-.1,length=0,lty=2)
mtext(expression(d),1,at=ncp,line=.5)

#-----------

# arrow for H1
arrows(cv+.1,.46,cv+2,.46,length=.15)
text(cv+.5,.48,'"H1"')

# arrow for H0
arrows(cv-.1,.46,cv-2,.46,length=.15)
text(cv-.5,.48,'"H0"')

#-----------

par(xpd=TRUE)

# indicate D0
arrows(D0,-.06,D0,-.005,length=.15,lty=1,lwd=1)
mtext(expression(D[0]),1,at=D0,line=3)

# indicate D1
arrows(D1,-.06,D1,-.005,length=.15,lty=1,lwd=1)
mtext(expression(D[1]),1,at=D1,line=3)

# indicate D
arrows(D,-.06,D,-.005,length=.15,lty=1,lwd=2)
mtext(expression(D),1,at=D,line=2.8)


#----------------------------------

dev.off()

