#!/usr/bin/env perl
#
# Modified from: https://github.com/lh3/psmc/blob/master/utils/ms2psmcfa.pl

use strict;
use warnings;
use Getopt::Std;

my %opts = (c=>0,l=>0,s=>100);
getopts('c:l:s:', \%opts);
my ($chrom, $len, $skip) = ($opts{c}, $opts{l}, $opts{s});
my @seq;

die "need chromosome name (-c param)" if ($chrom==0);
die "need chromosome length (-l param)" if ($len==0);

while (<>) {
  if (/^positions:/) {
	my @seq = ();
	my $l = int($len/$skip) + 1;
	$seq[$_] = 0 for (0 .. $l-1);

	chomp;
	my @segsites = split(/\s/);
	@segsites = @segsites[1..$#segsites];
	foreach (@segsites) {
		$seq[int($_ * $len / $skip)] = 1;
	}

	print ">$chrom";
	for my $i (0 .. $#seq) {
	  print "\n" if ($i % 60 == 0);
	  print ($seq[$i]?'K':'T');
	}
	print "\n";
	last;
  }
}
