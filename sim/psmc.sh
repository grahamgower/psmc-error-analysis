#!/bin/sh

source ./model.settings

p=0
for pparam in $pparams; do
	p=$((p+1))
	for m in 1 2 3; do
		$psmc \
			-p $pparam \
			-o $odir/m$m.p-x$p.psmc.out \
			$odir/m$m.psmc_input.psmcfa &
	done
	wait
done
