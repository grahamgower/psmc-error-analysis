#!/bin/sh

source ./model.settings
odir=./out-truncpsmc
idir=./out

for m in `seq $models`; do
	for x in $fraglevels; do
		mkdir -p $odir/m$m/x$x
		./frag_ms2smc-truncpsmc.pl \
			-x $x \
			-l $chrlen \
			-p $odir/m$m/x$x/ \
			< $idir/m$m/scrm.txt &
	done
done
wait

