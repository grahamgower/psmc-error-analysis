#!/bin/sh

source ./model.settings

p=0
for pparam in $pparams; do
	p=$((p+1))
	for m in 1 2 3; do
		$msmc \
			-p $pparam \
			--fixedRecombination \
			-t 4 \
			-o $odir/m$m.p-x$p.msmc1.out \
			$odir/m$m.chr*.msmc_input.txt &
	done
	wait
done
