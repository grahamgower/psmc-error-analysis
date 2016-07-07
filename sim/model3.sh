#!/bin/sh
#
# model 3: bottleneck with exponential population size increase

source ./model.settings
m=3

# Going backwards in time, we start with size N0, then at generation T1,
# start exponentially decreasing until we reach size N1 at generation T2,
# then the population instantaneously increases to size N2.
N0=1000000
N1=10000
N2=100000
T1=2000
T2=50000

t1=`perl -e "print $T1/(4*$N0)"`
t2=`perl -e "print $T2/(4*$N0)"`
alpha1=`perl -e "print -log($N1/$N0)/$t2"`
n2=`perl -e "print $N2/$N0"`

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
