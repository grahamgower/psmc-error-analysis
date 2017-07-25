#!/bin/sh

source ./model.settings

# msmc needs a boatload of files
ulimit -Sn `ulimit -Hn` || (echo "Cannot set 'ulimit -Sn'"; exit 1)

threads=1

for p in `seq $pops`; do
	for m in `seq $models`; do
		for s in `seq 0 $((samples-1))`; do
			for x in $fraglevels; do
				mkdir -p $odir/m$m/msmc1/p$p/s$s/x$x
				$msmc \
					-p $pparam \
					--verbose \
					-t $threads \
					-i 30 \
					-o $odir/m$m/msmc1/p$p/s$s/x$x/msmc1.out \
					$odir/m$m/p$p/x$x/msmc.sample$s.* \
					|| exit 1
			done
		done
	done &

done
wait
