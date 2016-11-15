#!/bin/sh

source ./model.settings
#odir=out-truncpsmc

for m in `seq $models`; do
	for s in `seq 0 $((samples-1))`; do
		for x in $fraglevels; do
			mkdir -p $odir/m$m/psmc/s$s/x$x
			$psmc \
				-p $pparam \
				-o $odir/m$m/psmc/s$s/x$x/psmc.out \
				$odir/m$m/x$x/psmc.sample$s
			done
	done &
done
wait
