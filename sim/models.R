
model_sim <- function(model, xlim=c(1e2, 1e6), step=100) {
	t <- seq(xlim[1], xlim[2], step)
	m <- data.frame(time=t, N=length(t), GOF=rep(0,length(t)))
	if (model == 1) {
		# model 1: constant population size
		m$N <- rep(10000, length(t))
	} else if (model == 2) {
		# model 2: exponential population size decrease
		N0 <- 1000
		N1 <- 10000
		T1 <- 2000
		T2 <- 5000
		alpha <- -log(N1/N0)/(T2-T1)
		m$N[1:(T1/step)] <- N0
		m$N[(T1/step+1):(T2/step)] <- N0*exp(-alpha * (t[(T1/step+1):(T2/step)]-T1))
		m$N[(T2/step+1):length(t)] <- N1
	} else if (model == 3) {
		# model 3: bottleneck with exponential population size increase
		N0 <- 100000
		N1 <- 1000
		N2 <- 100000
		T1 <- 2000
		T2 <- 50000
		alpha <- -log(N1/N0)/(T2-T1)
		m$N[1:(T1/step)] <- N0
		m$N[(T1/step+1):(T2/step)] <- N0*exp(-alpha * (t[(T1/step+1):(T2/step)]-T1))
		m$N[(T2/step+1):length(t)] <- N2
	}
	class(m) <- c("PSMC", class(m))
	return (m)
}

m1 <- model_sim(1)
m2 <- model_sim(2)
m3 <- model_sim(3)
