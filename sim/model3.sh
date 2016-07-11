#!/bin/sh
#
# model 3: bottleneck with exponential population size increase

source ./model.settings
m=3

# Going backwards in time, we start with size N0, then at generation T1,
# start exponentially decreasing until we reach size N1 at generation T2,
# then the population instantaneously increases to size N2.
N0=100000
N1=1000
N2=100000
T1=2000
T2=50000

t1=`perl -e "print $T1/(4*$N0)"`
t2=`perl -e "print $T2/(4*$N0)"`
alpha1=`perl -e "print -log($N1/$N0)/($t2-$t1)"`
n2=`perl -e "print $N2/$N0"`


mkdir -p $odir

chr=0
for size in $chrsizes; do
	chr=$((chr+1))

	theta=`perl -e "print $size*4*$N0*$mu"`
	recom=`perl -e "print $size*4*$N0*$mu/$mu_on_r"`

	$scrm \
		$samples \
		1 \
		-t $theta \
		-r $recom $size \
		-eG $t1 $alpha1 \
		-eN $t2 $n2 \
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
