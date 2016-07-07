#!/bin/sh
#
# model 2: exponential population size decrease

source ./model.settings
m=2

# Going backwards in time, we start with size N0, then at generation T1,
# start exponentially increasing until we reach size N1 at generation T2.
N0=10000
N1=1000000
T1=2000
T2=50000

t1=`perl -e "print $T1/(4*$N0)"`
t2=`perl -e "print $T2/(4*$N0)"`
alpha1=`perl -e "print -log($N1/$N0)/($t2-$t1)"`
n2=`perl -e "print $N1/$N0"`

theta=`perl -e "print 4*$N0*$mu"`
recom=`perl -e "print 4*$N0*$mu/$mu_on_r"`

chr=0
for size in $chrsizes; do
	chr=$((chr+1))

	$macs \
		$samples \
		$size \
		-t $theta \
		-r $recom \
		-h 1e2 \
		-eG $t1 $alpha1 \
		-eG $t2 0 \
		-eN $t2 $n2 \
		2> m$m.chr$chr.trees.txt \
		1> m$m.chr$chr.sites.txt
done
