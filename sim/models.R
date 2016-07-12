#!/usr/bin/env Rscript

# time, in generations
t = 1:2e5

# model 1: constant population size
m1.N = rep(10000, length(t))

# model 2: exponential population size decrease
m2.N0 = 1000
m2.N1 = 10000
m2.T1 = 2000
m2.T2 = 5000
m2.alpha1 = -log(m2.N1/m2.N0)/(m2.T2-m2.T1)
m2.N = vector(length=length(t))
m2.N[1:m2.T1] = m2.N0
m2.N[m2.T1:m2.T2] = m2.N0*exp(-m2.alpha1 * (t[m2.T1:m2.T2]-m2.T1))
m2.N[m2.T2:length(t)] = m2.N1

# model 3: bottleneck with exponential population size increase
m3.N0 = 100000
m3.N1 = 1000
m3.N2 = 100000
m3.T1 = 2000
m3.T2 = 50000
m3.alpha1 = -log(m3.N1/m3.N0)/(m3.T2-m3.T1)
m3.N = vector(length=length(t))
m3.N[1:m3.T1] = m3.N0
m3.N[m3.T1:m3.T2] = m3.N0*exp(-m3.alpha1 * (t[m3.T1:m3.T2]-m3.T1))
m3.N[m3.T2:length(t)] = m3.N2

vmin = function(v) v[which.min(v)]
vmax = function(v) v[which.max(v)]
ymin = vmin(c(m1.N, m2.N0, m2.N1, m3.N0, m3.N1, m3.N2))
ymax = vmax(c(m1.N, m2.N0, m2.N1, m3.N0, m3.N1, m3.N2))
xmin = 100
xmax = 1e5


llog="xy"
lwd = 2
lty = c(1,2,3)
col = c("red", "blue", "black")
pdf("models.pdf", width=12, height=8)
plot(t, m1.N, type='l', lty=lty[1], lwd=lwd, col=col[1], log=llog,
     ylim=c(ymin, ymax), xlim=c(xmin, xmax), xaxt="n", yaxt="n", xlab="", ylab="")
par(new=T)
plot(t, m2.N, type='l', lty=lty[2], lwd=lwd, col=col[2], log=llog,
     ylim=c(ymin, ymax), xlim=c(xmin, xmax), xaxt="n", yaxt="n", xlab="", ylab="")
par(new=T)
plot(t, m3.N, type='l', lty=lty[3], lwd=lwd, col=col[3], log=llog,
     ylim=c(ymin, ymax), xlim=c(xmin, xmax), xaxt="n", yaxt="n", xlab="", ylab="")
xticks=c(seq(1e2,1e3,1e2), seq(1e3,1e4,1e3), seq(1e4,1e5,1e4))
yticks=c(seq(1e3,1e4,1e3), seq(1e4,1e5,1e4), seq(1e5,1e6,1e5))
axis(1, xticks, labels=F)
axis(2, yticks, labels=F)
axis(4, yticks, labels=F)
axis(1, c(1e2,1e3,1e4,1e5), lwd.ticks=3)
axis(2, c(1e3,1e4,1e5,1e6), lwd.ticks=3)
axis(4, c(1e3,1e4,1e5,1e6), lwd.ticks=3, labels=F)
title("Demographic models for PSMC/PSMC' error analysis",
      xlab="time (generations)", ylab="Ne")
legend("topleft", c("model 1", "model 2", "model 3"), lwd=2, lty=lty, col=col)
invisible(dev.off())
