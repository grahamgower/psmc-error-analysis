#!/usr/bin/env perl
# Fragment ms simulated data, using an empirical scaffold length distribution,
# outputting psmc and msmc input files.

use strict;
use warnings;
use Getopt::Std;
use Statistics::Basic 'median';

my %opts = (c=>"", f=>"", l=>1e6, p=>"out");
getopts("c:f:l:", \%opts);
my ($chr, $faidx, $chrlen, $prefix) = ($opts{c}, $opts{f}, int($opts{l}), $opts{p});
my (@fraglens, $median_fraglen);
my @segsites;

my $skip = 100; # psmc bin size
my $minlen = 1e6;
srand(31415);

if ($faidx eq "" || $chr eq "" || $chrlen < 1) {
    print STDERR "Usage: $0 -f frag.fai -c chr -l chrlen -p outprefix ms.txt\n";
    exit 1;
}

open(my $f, $faidx) or die "Err: cannot open $faidx: $!";
while (<$f>) {
    chomp;
    my @fields = split(/\s/);
    my $len = int($fields[1]);
    if ($len < $minlen) {
        next
    }
    if ($len > $chrlen) {
        die "$faidx: line $.: scaffold length ($len) > chrlen ($chrlen)";
    }
    push(@fraglens, $len);
}
close($f);

$median_fraglen = median(@fraglens);

while (<>) {
  if (/^positions:/) {
	chomp;
	my @fields = split(/\s/);
	@segsites = map { int($_*$chrlen) } @fields[1..$#fields];
        last;
    }
}

my ($from, $to) = (0, 0);
my $i = 0;

open(my $pfh, ">", "$prefix.psmc")
    or die "cannot open $prefix.psmc for writing: $!";

while ($to+$median_fraglen < $chrlen) {
    my $len;
    do {
        $len = $fraglens[rand(@fraglens)];
    } while ($to+$len > $chrlen);

    $from = $to+1;
    $to += $len;
    $i++;
    #print("${chr}_${i}\t$from\t$to\n");

    my $k = 0;
    while ($segsites[$k] <= $to) {
        $k++;
    }
    my @sites = splice(@segsites, 0, $k);

    # psmc
    my @seq = ();
    $seq[$_] = 0 for (0 .. int($len/$skip));
    $seq[int(($_-$from)/$skip)] = 1 for (@sites);
    print $pfh ">${chr}_$i";
    for my $j (0 .. $#seq) {
        print ($pfh "\n") if ($j % 60 == 0);
        print ($pfh $seq[$j]?'K':'T');
    }
    print $pfh "\n";

    # msmc
    open(my $mfh, ">", "$prefix.msmc.${chr}_$i")
        or die "cannot open $prefix.msmc.${chr}_$i for writing: $!";
    my $last = 1;
    for (@sites) {
        my $gap = $_ - $last +1;
        print $mfh "${chr}_$i\t$_\t$gap\t10,01\n";
        $last = $_;
    }
    close($mfh);
}

close($pfh);
