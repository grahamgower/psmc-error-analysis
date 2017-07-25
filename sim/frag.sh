#!/bin/sh

source ./model.settings

do_frag()
{
	odir=$1
	chrlen=$2
	m=$3
	p=$4
	x=$5
	mkdir -p $odir/m$m/p$p/x$x
	./frag_ms2smc.pl \
		-x $x \
		-l $chrlen \
		-p $odir/m$m/p$p/x$x/ \
		< $odir/m$m/p$p/scrm.txt
}
export -f do_frag

for p in `seq $pops`; do
	for m in `seq $models`; do
		for x in $fraglevels; do
			echo do_frag $odir $chrlen $m $p $x
		done
	done
done | parallel -j 8
