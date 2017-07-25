#!/bin/sh

source ./model.settings

for p in `seq $pops`; do
	for m in `seq $models`; do
		for s in `seq 0 $((samples-1))`; do
			for x in $fraglevels; do
				mkdir -p $odir/m$m/psmc/p$p/s$s/x$x
				mkdir -p $odir/m$m/psmc-trunc/p$p/s$s/x$x
				echo $psmc \
					-p $pparam \
					-o $odir/m$m/psmc/p$p/s$s/x$x/psmc.out \
					$odir/m$m/p$p/x$x/psmc.sample$s
				echo $psmc \
					-p $pparam \
					-o $odir/m$m/psmc-trunc/p$p/s$s/x$x/psmc.out \
					$odir/m$m/p$p/x$x/psmc-trunc.sample$s
			done
		done
	done
done | parallel -j 8
