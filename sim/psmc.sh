#!/bin/sh

source ./model.settings

pparams="
4*1+25*2+4*1+6*1
1*2+15*1+1*2
10*1+15*2
"

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
