#!/bin/sh

source ./model.settings

# number of haploid samples
hsamples=`perl -e "print 2*$samples"`

# model 1: constant population size
model1() {
	m=1
	N=10000

	mkdir -p $odir/m$m

	theta=`perl -e "print $chrlen*4*$N*$mu"`
	recom=`perl -e "print $chrlen*4*$N*$mu/$mu_on_r"`

	$scrm \
		$hsamples \
		1 \
		-t $theta \
		-r $recom $chrlen \
		> $odir/m$m/scrm.txt
}

# model 2: exponential population size decrease
model2() {
	m=2

	# Going backwards in time, we start with size N0, then at generation T1,
	# start exponentially increasing until we reach size N1 at generation T2.
	N0=1000
	N1=10000
	T1=2000
	T2=5000

	t1=`perl -e "print $T1/(4*$N0)"`
	t2=`perl -e "print $T2/(4*$N0)"`
	alpha1=`perl -e "print -log($N1/$N0)/($t2-$t1)"`
	n2=`perl -e "print $N1/$N0"`

	theta=`perl -e "print $chrlen*4*$N0*$mu"`
	recom=`perl -e "print $chrlen*4*$N0*$mu/$mu_on_r"`

	mkdir -p $odir/m$m

	$scrm \
		$hsamples \
		1 \
		-t $theta \
		-r $recom $chrlen \
		-eG $t1 $alpha1 \
		-eG $t2 0.0 \
		> $odir/m$m/scrm.txt
}


# model 3: bottleneck with exponential population size increase
model3() {
	m=3

	# Going backwards in time, we start with size N0, then at generation T1,
	# start exponentially decreasing until we reach size N1 at generation T2,
	# then the population instantaneously increases to size N2.
	N0=100000
	N1=1000
	N2=100000
	T1=2000
	T2=10000

	t1=`perl -e "print $T1/(4*$N0)"`
	t2=`perl -e "print $T2/(4*$N0)"`
	alpha1=`perl -e "print -log($N1/$N0)/($t2-$t1)"`
	n2=`perl -e "print $N2/$N0"`

	theta=`perl -e "print $chrlen*4*$N0*$mu"`
	recom=`perl -e "print $chrlen*4*$N0*$mu/$mu_on_r"`

	mkdir -p $odir/m$m

	$scrm \
		$hsamples \
		1 \
		-t $theta \
		-r $recom $chrlen \
		-eG $t1 $alpha1 \
		-eN $t2 $n2 \
		> $odir/m$m/scrm.txt
}

model1 &
model2 &
model3 &
wait
