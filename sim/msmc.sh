#!/bin/sh

source ./model.settings

# msmc needs a boatload of files
ulimit -Sn 16384 || (echo "Cannot set 'ulimit -Sn 16384'"; exit 1)

threads=4

for m in `seq $models`; do
	for s in `seq 0 $((samples-1))`; do
		for x in $fraglevels; do
			mkdir -p $odir/m$m/msmc1/s$s/x$x
			$msmc \
				-p $pparam \
				--verbose \
				-t $threads \
				-i 30 \
				-o $odir/m$m/msmc1/s$s/x$x/msmc1.out \
				$odir/m$m/x$x/msmc.sample$s.*
		done
	done &
done
wait
