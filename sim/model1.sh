#!/bin/sh
#
# model 1: constant population size

source ./model.settings
m=1

N=10000
theta=`perl -e "print 4*$N*$mu"`
recom=`perl -e "print 4*$N*$mu/$mu_on_r"`

chr=0
for size in $chrsizes; do
	chr=$((chr+1))

	$macs \
		$samples \
		$size \
		-t $theta \
		-r $recom \
		-h 1e2 \
		2> m$m.chr$chr.trees.txt \
		1> m$m.chr$chr.sites.txt
done
