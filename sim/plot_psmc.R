#!/usr/bin/env Rscript

printf <- function(...) cat(sprintf(...))

# Parse PSMC output, scaling the for mutation rate mu and generations g,
#  mu has units substitutions/base/gen.
# Returns time and population size for the specified round.
# n is the number of time bins PSMC used, which are further grouped
#  with the -p param to PSMC.
# s is the number of nucleotides encoded as one char in the PSMC input.
psmc.input = function(file, mu, g, round=20, n=64, s=100) {
	d = list() # list of data frames, one for each round
	f = file(file, "r")
	lineno = 0

	while (length(line <- readLines(f, n=1)) > 0) {
		lineno = lineno + 1
		if (grepl("^CC", line) || grepl("^MM", line))
			next
		fields = strsplit(line, "\t")[[1]]
		if (grepl("^RD", line)) {
			rd = as.integer(fields[2])+1
			d[[rd]] = data.frame(t=vector("numeric",length=n), lambda=vector("numeric",length=n))
				#pi = vector("numeric",length=n),
				#sum.Akl = vector("numeric",length=n),
				#Akk = vector("numeric",length=n))
		#} else if (grepl("^LK", line)) {
		#} else if (grepl("^QD", line)) {
		#} else if (grepl("^RI", line)) {
		} else if (grepl("^TR", line)) {
			theta = as.double(fields[2])
			rho = as.double(fields[3])
			d[[rd]]$theta = theta
			d[[rd]]$rho = rho
		#} else if (grepl("^MT", line)) {
		} else if (grepl("^RS", line)) {
			k = as.integer(fields[2])+1
			if (k > n) {
				printf("%s: line %d: n=%d is too small.\n",
				       file, lineno, n)
				break
			}
			d[[rd]]$t[k] = as.double(fields[3])
			d[[rd]]$lambda[k] = as.double(fields[4])
			#d[[rd]]$pi[k] = as.double(fields[5])
			#d[[rd]]$sum.Akl[k] = as.double(fields[6])
			#d[[rd]]$Akk[k] = as.double(fields[7])
		#} else if (grepl("^PA", line)) {
		}
	}

	close(f)

	m = data.frame(t=d[[round]]$t, lambda=d[[round]]$lambda)
	N0 = d[[round]]$theta[1] / (4*mu) / s
	m$time = 2*N0 * m$t / g
	m$N = N0 * m$lambda
	return (m)
}

# Parse MSMC final.txt, scaling the for mutation rate mu and generations g,
#  mu has units substitutions/base/gen.
msmc.input = function(file, mu, g) {
	m = read.table(file, header=T)
	m$time = g * m$left_time_boundary / mu
	m$N = 1.0 / (2*mu*m$lambda)
	return (m)
}

# find min/max of a vector - I can't believe R doesn't have this
vmin = function(v) v[which.min(v)]
vmax = function(v) v[which.max(v)]

# simulated model
model.sim = function(model, step=100) {
	t = seq(1,1e6,step)
	m = data.frame(time=t, N=length(t))
	if (model == 1) {
		# model 1: constant population size
		m$N = rep(10000, length(t))
	} else if (model == 2) {
		# model 2: exponential population size decrease
		N0 = 1000
		N1 = 10000
		T1 = 2000
		T2 = 5000
		alpha = -log(N1/N0)/(T2-T1)
		m$N[1:(T1/step)] = N0
		m$N[(T1/step+1):(T2/step)] = N0*exp(-alpha * (t[(T1/step+1):(T2/step)]-T1))
		m$N[(T2/step+1):length(t)] = N1
	} else if (model == 3) {
		# model 3: bottleneck with exponential population size increase
		N0 = 100000
		N1 = 1000
		N2 = 100000
		T1 = 2000
		T2 = 50000
		alpha = -log(N1/N0)/(T2-T1)
		m$N[1:(T1/step)] = N0
		m$N[(T1/step+1):(T2/step)] = N0*exp(-alpha * (t[(T1/step+1):(T2/step)]-T1))
		m$N[(T2/step+1):length(t)] = N2
	}
	return (m)
}

# mutations per base per generation
mu=1.25e-8

args = commandArgs(trailingOnly=T)
if (length(args) != 5) {
	arg0 = strsplit(commandArgs()[4],"=")[[1]][2]
	printf("usage: %s model psmc.out msmc.out out.pdf title\n", arg0)
	quit("no", 1)
}

d1 = model.sim(as.integer(args[1]))
d2 = psmc.input(args[2], mu, 1)
d3 = msmc.input(args[3], mu, 1)

# bounds for the plot
xmin = 1e2 #vmin(d1$time[2:length(d1$time)])
xmax = 5e5 #vmax(d1$time)
ymin = 1e3 #vmin(d1$N)
ymax = 1e6 #vmax(d1$N[1:(length(d1$N)-1)])

pdf(args[4], width=12, height=8)
lw=3
lty = c(3,1,2)
col = c("black", "red", "blue")
plot(d1$time, d1$N, log="xy", ylim=c(ymin,ymax), xlim=c(xmin,xmax),
     type="n", xlab="", ylab="", xaxt="n", yaxt="n")
lines(d1$time, d1$N, type="s", lty=lty[1], col=col[1], lw=lw)
lines(d2$time, d2$N, type="s", lty=lty[2], col=col[2], lw=lw)
lines(d3$time, d3$N, type="s", lty=lty[3], col=col[3], lw=lw)

legend("topleft", c("simulation", "PSMC", "MSMC"), lwd=lw, lty=lty, col=col)
title(args[5], xlab="time (generations in the past)", ylab="Ne")

# sensible tick marks and labels for log-log axes
ticks=c(seq(1e2,1e3,1e2), seq(1e3,1e4,1e3), seq(1e4,1e5,1e4), seq(1e5,1e6,1e5))
bdticks=c(1e2,1e3,1e4,1e5,1e6)
axis(1, ticks, labels=F)
axis(2, ticks, labels=F)
axis(4, ticks, labels=F)
axis(1, bdticks, lwd.ticks=3)
axis(2, bdticks, lwd.ticks=3)
axis(4, bdticks, lwd.ticks=3, labels=F)

invisible(dev.off())
