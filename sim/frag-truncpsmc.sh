#!/bin/sh

odir=./out-truncpsmc
idir=./out

chrlen=100000000

# levels of fragmentation
xx="100000000 10000000 1000000 100000 10000"

for m in 1 2 3; do
	for x in $xx; do
		mkdir -p $odir/m$m/x$x
		./frag_ms2smc-truncpsmc.pl \
			-x $x \
			-l $chrlen \
			-p $odir/m$m/x$x/ \
			< $idir/m$m/scrm.txt &
	done
done
wait
