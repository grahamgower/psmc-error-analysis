#!/bin/sh
#
# We use lots of small files, so try a ramdisk instead. `zram' provides
# compression on RAM based block devices, allowing this to fit.
# https://www.kernel.org/doc/Documentation/blockdev/zram.txt

devno=0
dev=/dev/zram$devno
mnt=/mnt/zram
user=grg
group=acad
ddir=/sys/block/zram$devno/
disksize=48G
memlimit=32G

modprobe zram
umount $mnt

echo 1 > $ddir/reset

echo $disksize > $ddir/disksize
echo $memlimit > $ddir/mem_limit

mkfs.ext4 $dev
mkdir -p $mnt
mount -o rw,noatime,nodev,nosuid $dev $mnt

chown ${user}:${group} $mnt
