#!/bin/sh

source ./model.settings

for m in `seq $models`; do
	for x in $fraglevels; do
		mkdir -p $odir/m$m/x$x
		./frag_ms2smc.pl \
			-x $x \
			-l $chrlen \
			-p $odir/m$m/x$x/ \
			< $odir/m$m/scrm.txt &
	done
done
wait

