#!/bin/sh
#
# model 1: constant population size

source ./model.settings

m=1
N=10000

mkdir -p $odir

chr=0
for size in $chrsizes; do
	chr=$((chr+1))

	theta=`perl -e "print $size*4*$N*$mu"`
	recom=`perl -e "print $size*4*$N*$mu/$mu_on_r"`

	$scrm \
		$samples \
		1 \
		-t $theta \
		-r $recom $size \
		> $odir/m$m.chr$chr.scrm.txt

	$ms2multihetsep $chr $size \
		< $odir/m$m.chr$chr.scrm.txt \
		> $odir/m$m.chr$chr.msmc_input.txt
done

chr=0
for size in $chrsizes; do
	chr=$((chr+1))
	./ms2psmcfa.pl -c $chr -l $size \
		< $odir/m$m.chr$chr.scrm.txt
done > $odir/m$m.psmc_input.psmcfa
