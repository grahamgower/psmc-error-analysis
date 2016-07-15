#!/usr/bin/env Rscript

pdf("scafflens_qqexp.pdf")

for (file in Sys.glob("faidx/*.fai")) {
	fai = read.table(file)
	y = fai[fai$V2>1e6,]$V2
	if (length(y) == 0) {
		next
	}
	qqplot(qexp(ppoints(length(y))), y)
	qqline(y, distribution=qexp)
	mtext(file)
}

invisible(dev.off())
