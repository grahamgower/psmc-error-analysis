#!/bin/sh

source ./model.settings

pdflist=""
m=0
for model in "constant population size" "exponential decrease" "bottleneck with recovery"; do
	m=$((m+1))
	p=0
	for pparam in $pparams; do
		p=$((p+1))
		./plot_psmc.R \
			$m \
			out/m$m.p-x$p.psmc.out \
			out/m$m.p-x$p.msmc1.out.final.txt \
			m$m.p-x$p.pdf \
			"model $m: $model: -p $pparam"
		pdflist="$pdflist m$m.p-x$p.pdf"
	done
done

#./models.R
pdflist="models.pdf $pdflist"

pdfunite $pdflist psmc-models.pdf
