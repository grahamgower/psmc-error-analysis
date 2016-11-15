
source('psmcr/R/psmcr.R')

# mutations per base per generation
mu <- 1.25e-8
# generation time
g <- 1

models <- 1:3
fragsizes <- c("100000000", "10000000", "1000000", "100000", "10000")
nsamples <- 100

d = data.frame(prog=factor(levels=c("PSMC", "PSMC-trunc","MSMC1")),
		model=factor(levels=models),
		fragsize=vector("numeric"),
		sample=vector("numeric"),
		rd=vector("numeric"), # EM round
		GOF=vector("numeric"),
		ll=vector("numeric"),
		theta=vector("numeric"), # PSMC only, EM inferred
		rho=vector("numeric"), # EM inferred
		T_max=vector("numeric"), # PSMC only, EM inferred
		k=vector("numeric"), # time bin
		t=vector("numeric"), # time in 2*N0 generations
		lambda=vector("numeric"),
		pi <- vector("numeric"), # PSMC only, actually \pi_k*C_pi
		sigma <- vector("numeric"),
		post_sigma <- vector("numeric"),
		time <- vector("numeric"), # time in years
		N <- vector("numeric")) # effective pop. size

for (m in models) {
	for (x in fragsizes) {
		for (s in 1:nsamples) {
			psmc_file <- sprintf("out/m%d/psmc/s%d/x%s/psmc.out", m, s-1, x)
			psmc <- lastround(parse_psmc(psmc_file, mu, g))
			d = rbind(d, cbind(prog="PSMC", model=m, fragsize=as.integer(x), sample=s, psmc))

			psmc_file <- sprintf("out-truncpsmc/m%d/psmc/s%d/x%s/psmc.out", m, s-1, x)
			psmc <- lastround(parse_psmc(psmc_file, mu, g))
			d = rbind(d, cbind(prog="PSMC-trunc", model=m, fragsize=as.integer(x), sample=s, psmc))

			msmc1_filebase <- sprintf("out/m%d/msmc1/s%d/x%s/msmc1.out", m, s-1, x)
			msmc1 = lastround(parse_msmc(msmc1_filebase, mu, g))
			d = rbind(d, cbind(prog="MSMC1", model=m, fragsize=as.integer(x), sample=s, msmc1))
		}
	}
}

fn <- sprintf("psmc-error-analysis_%s.tsv.gz", strftime(Sys.time(), "%Y%m%d"))
gz1 <- gzfile(fn, "w")
write.table(d, file=gz1, sep="\t", quote=F, row.names=F)
close(gz1)
