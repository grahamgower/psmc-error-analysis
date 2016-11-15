# Parse a pattern like "4+5*3+4".
# Stolen from psmc/cli.c:psmc_parse_pattern()
psmc_parse_pattern <- function(free_lambdas, pattern) {
	lambdas <- list()
	x <- 1
	y <- 1
	k <- 1
	p <- q <- 1
	while (p <= nchar(pattern)) {
		a <- substr(pattern, p, p)
		if (a == '+') {
			l <- as.integer(substr(pattern, q, p-1))
			p <- p + 1
			q <- p
			for (i in 1:k) {
				for (j in 1:l) {
					lambdas[x] <- free_lambdas[y]
					x <- x + 1
				}
				y <- y + 1
			}
			k <- 1
		} else if (a == '*') {
			k <- as.integer(substr(pattern, q, p-1)) # number of repeats
			p <- p + 1
			q <- p
		} else {
			p <- p + 1
		}
	}

	if (q < p) {
		l <- as.integer(substr(pattern, q, p-1))
		for (i in 1:k) {
			for (j in 1:l) {
				lambdas[x] <- free_lambdas[y]
				x <- x + 1
			}
			y <- y + 1
		}
	}

	return (as.vector(lambdas, "numeric"))
}

# Parse PSMC output, scaling for mutation rate mu and generations g,
#  mu has units substitutions/base/generation.
# s is the number of nucleotides encoded as one char in the PSMC input.
parse_psmc <- function(file, mu, g, s=100) {
	d <- data.frame(GOF=vector("numeric"),
			ll=vector("numeric"),
			theta=vector("numeric"),
			rho=vector("numeric"),
			T_max=vector("numeric"))
	dlist <- list()

	diverge <- 0
	f <- file(file, "r")
	lineno <- 0

	while (length(line <- readLines(f, n=1)) > 0) {
		lineno <- lineno + 1
		fields <- strsplit(line, "\t")[[1]]
		if (grepl("^RD", line)) {
			rd <- as.integer(fields[2])+1
			d[rd,] <- rep(NA, length(d))
			dlist[[rd]] <- data.frame(k=vector("numeric"),
						t=vector("numeric"),
						lambda=vector("numeric"),
						pi=vector("numeric"),
						sigma=vector("numeric"),
						post_sigma=vector("numeric"))
		} else if (grepl("^MM.*C_pi", line)){
			tmp1 <- strsplit(fields[2], " ")[[1]][2]
			tmp2 <- substr(tmp1, start=1, stop=nchar(tmp1)-1)
			d$C_pi[rd] <- as.double(tmp2)
		} else if (grepl("^LK", line)) {
			d$ll[rd] <- as.double(fields[2])
		#} else if (grepl("^QD", line)) {
		} else if (grepl("^RI", line)) {
			d$GOF[rd] <- as.double(fields[2])
		} else if (grepl("^TR", line)) {
			#d$theta[rd] <- as.double(fields[2])
			#d$rho[rd] <- as.double(fields[3])
		} else if (grepl("^MT", line)) {
			#d$T_max[rd] <- as.double(fields[2])
		} else if (grepl("^DT", line)) {
			diverge <- 1
			#d$dt[rd] <- as.double(fields[2])
		} else if (grepl("^RS", line)) {
			k <- as.integer(fields[2])
			dlist[[rd]][k+1,] <- lapply(fields[2:7], as.double)
		} else if (grepl("^PA", line)) {
			# high precision: %.9lf vs. %.6lf for other lines
			# pparam, theta, rho, T_max, {lambda_k for each k}
			fields2 <- strsplit(fields[2], " ")[[1]]
			pattern <- fields2[1]
			d$theta[rd] <- as.double(fields2[2])
			d$rho[rd] <- as.double(fields2[3])
			d$T_max[rd] <- as.double(fields2[4])

			free_lambdas <- lapply(fields2[5:(length(fields2)-diverge)], as.double)
			dlist[[rd]]$lambda <- psmc_parse_pattern(free_lambdas, pattern)
		}
	}

	close(f)

	d$N0 <- d$theta / (4*mu) / s
	d$C_sigma <- 1.0 / (d$C_pi * d$rho) + 0.5

	psmc <- data.frame(
			rd=vector("numeric"),
			GOF=vector("numeric"),
			ll=vector("numeric"),
			theta=vector("numeric"),
			rho=vector("numeric"),
			T_max=vector("numeric"),
			k=vector("numeric"),
			t=vector("numeric"),
			lambda=vector("numeric"),
			pi=vector("numeric"),
			sigma=vector("numeric"),
			post_sigma=vector("numeric"),
			time=vector("numeric"),
			N=vector("numeric"))

	i <- 0
	for (r in 1:rd) {
		dlist[[r]]$time <- 2*d$N0[r] * dlist[[r]]$t / g
		dlist[[r]]$N <- d$N0[r] * dlist[[r]]$lambda

		# RI line above has greater precision
		#d$GOF[r] <- sum(dlist[[r]]$sigma * log(dlist[[r]]$sigma/dlist[[r]]$post_sigma), na.rm=T)

		for (kk in 1:(k+1)) {
			i <- i+1
			psmc[i,] <- c(r-1, d$GOF[r], d$ll[r], d$theta[r], d$rho[r], d$T_max[r],
				      kk-1, dlist[[r]]$t[kk], dlist[[r]]$lambda[kk], dlist[[r]]$pi[kk], dlist[[r]]$sigma[kk], dlist[[r]]$post_sigma[kk], dlist[[r]]$time[kk], dlist[[r]]$N[kk])
		}
	}

	class(psmc) <- c("PSMC", class(psmc))
	return (psmc)
}

# Parse MSMC output, with output prefix pfx, scaling for mutation rate mu
# and generations g. mu has units substitutions/base/gen.
parse_msmc <- function(pfx, mu, g) {

	#file.final <- sprintf("%s.final.txt", pfx)
	#f <- read.table(file.final, header=T)
	#f$time <- g * f$left_time_boundary / mu
	#f$N <- 1.0 / (2*mu*f$lambda)

	Vlist <- list()
	#Mlist <- list()
	i <- 0
	file.loop_EM <- sprintf("%s.loop_%d.expectationMatrix.txt", pfx, i);
	while (file.exists(file.loop_EM)) {
		i <- i+1
		# posterior expectation of state, from HMM
		Vlist[[i]] <- as.vector(read.table(file.loop_EM, header=F, nrows=1))
		# posterior expectation of state transition, from HMM
		#Mlist[[i]] <- as.matrix(read.table(file.loop_EM, header=F, skip=1))
		file.loop_EM <- sprintf("%s.loop_%d.expectationMatrix.txt", pfx, i);
	}

	file.loop <- sprintf("%s.loop.txt", pfx)
	loop <- read.table(file.loop, header=F, stringsAsFactors=F,
			  col.names=c("rho", "ll", "time_boundaries", "lambdas"))
	time_boundaries <- lapply(strsplit(loop$time_boundaries, ","), as.double)
	lambdas <- lapply(strsplit(loop$lambdas, ","), as.double)
	rounds <- length(lambdas)

	if (i > 0 && i != rounds) {
		stop(sprintf("%s has %d rounds, but %d *.loop_*.expectationMatrix.txt files found\n",
			     file.loop, rounds, i))
	}

	n <- length(lambdas[[1]])
	d <- data.frame(rho=loop$rho, ll=loop$ll)
	dlist <- list()

	msmc <- data.frame(
			rd=vector("numeric"),
			GOF=vector("numeric"),
			ll=vector("numeric"),
			theta=vector("numeric"),
			rho=vector("numeric"),
			T_max=vector("numeric"),
			k=vector("numeric"),
			t=vector("numeric"),
			lambda=vector("numeric"),
			pi=vector("numeric"),
			sigma=vector("numeric"),
			post_sigma=vector("numeric"),
			time=vector("numeric"),
			N=vector("numeric"))

	j <- 0
	for (r in 1:rounds) {
		dlist[[r]] <- data.frame(t=time_boundaries[[r]][1:n],
					lambda=lambdas[[r]][1:n])

		time_boundaries[[r]][n+1] <- 1000 # clobber Inf
		tau <- diff(time_boundaries[[r]])
		sum_area <- vector("numeric", n)
		sum_area[1] <- 0
		for (k in 2:n) {
			sum_area[k] <- sum_area[k-1] + tau[k-1]*dlist[[r]]$lambda[k-1]
		}
		sigma_ <- dlist[[r]]$lambda * exp(-sum_area)
		dlist[[r]]$sigma <- sigma_ / sum(sigma_)

		if (i > 0) {
			dlist[[r]]$post_sigma <- t(Vlist[[r]] / sum(Vlist[[r]]))
			d$GOF[r] <- sum(dlist[[r]]$sigma * log(dlist[[r]]$sigma / dlist[[r]]$post_sigma), na.rm=T)
		} else {
			dlist[[r]]$post_sigma <- rep(NA, n)
			d$GOF[r] <- NA
		}

		dlist[[r]]$time <- g * dlist[[r]]$t / mu
		dlist[[r]]$N <- 1.0 / (2*mu*dlist[[r]]$lambda)

		for (k in 1:n) {
			j <- j+1
			msmc[j,] <- c(r, d$GOF[r], d$ll[r], NA, d$rho[r], NA,
				      k-1, dlist[[r]]$t[k], dlist[[r]]$lambda[k], NA, dlist[[r]]$sigma[k], dlist[[r]]$post_sigma[k], dlist[[r]]$time[k], dlist[[r]]$N[k])
		}
	}

	class(msmc) <- c("PSMC", class(msmc))
	return (msmc)
}

lastround <- function(psmc) psmc[psmc$rd==max(psmc$rd),]

bestfit <- function(psmc) {
	if ("GOF" %in% colnames(psmc)) {
		return (psmc[psmc$GOF==min(psmc$GOF),])
	} else {
		stop("GOF undefined, try running msmc with --verbose")
		# use the last round for MSMC runs without --verbose
		#return (lastround(psmc))
	}
}

plot.PSMC <- function(psmc, ...) {
	b <- bestfit(psmc)
	plot(b$time, b$N, type='s', log='xy', ...)
}

ggplot.PSMC <- function(psmc, mapping=aes(time,N), ...) {
	require(ggplot2)
	b <- bestfit(psmc)
	mbreaks <- function(a) {x<-c(); y<-10; while(y<a[1]) y<-y*10; while(y<=a[2]) {x<-c(x,seq(y*2,y*10,y)); y<-y*10} ; return(x)}
	lfunc <- function(a) parse(text=sprintf("10^%.0f",log10(a)))
	m <- NextMethod(data=b, mapping=mapping, ...) +
		geom_step() +
		scale_x_log10(minor_breaks=mbreaks, labels=lfunc) +
		scale_y_log10(minor_breaks=mbreaks, labels=lfunc) +
		annotation_logticks()
	return (m)
}
