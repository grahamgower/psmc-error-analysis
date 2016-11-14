#!/usr/bin/env perl
# Fragment ms simulated data, using an empirical scaffold length distribution,
# outputting psmc and msmc input files.

use strict;
use warnings;
use Getopt::Std;

my %opts = (l=>-1, p=>"out", x=>1);
getopts("l:p:x:", \%opts);
my ($chrlen, $prefix, $fraglevel) = (int($opts{l}), $opts{p}, $opts{x});
my (@segsites, @haplotypes);

my $skip = 100; # psmc bin size
my $minlen = 10*$skip;
srand(31415*$fraglevel);

if ($fraglevel < $minlen) {
    print STDERR "Error: -x $fraglevel too low, must be >= $minlen";
    exit 1;
}

if ($chrlen < $fraglevel) {
    print STDERR "Usage: $0 -l chrlen -x fraglevel -p outprefix ms.txt\n";
    exit 1;
}

while (<>) {
    if (/^positions:/) {
	chomp;
	my @fields = split(/\s/);
	@segsites = map { int($_*$chrlen) } @fields[1..$#fields];
    } elsif (@segsites) {
        if (/^$/) {
            last;
        }
        chomp;
        push(@haplotypes, $_);
    }
}

for my $sample (0 .. $#haplotypes/2) {

    my ($from, $to) = (0, 0);
    my $i = 0;

    my @ss = @segsites;
    my $hap1 = $haplotypes[2*$sample];
    my $hap2 = $haplotypes[2*$sample+1];

    open(my $pfh, ">", "$prefix/psmc.sample$sample")
        or die "cannot open $prefix/psmc.sample$sample for writing: $!";

    while ($to < $chrlen) {
        my $len = $fraglevel;
        if ($to+$len > $chrlen) {
            # just use remainder of the chromosome
            $len = $chrlen - $to;
        }

        $from = $to+1;
        $to += $len;
        $i++;
        #print("${i}\t$from\t$to\n");

        my $k = 0;
        while ($k < @ss and $ss[$k] <= $to) {
            $k++;
        }

        my @sites = splice(@ss, 0, $k);
        my $h1 = substr($hap1, 0, $k, "");
        my $h2 = substr($hap2, 0, $k, "");

        # psmc
        my $noseg = 1;
        my $lastpos;
        my @seq = ();
        $seq[$_] = 0 for (0 .. int($len/$skip));
        for my $x (0 .. $#sites) {
            if (substr($h1,$x,1) != substr($h2,$x,1)) {
                $seq[int(($sites[$x]-$from)/$skip)] = 1;
                $lastpos = int(($sites[$x]-$from)/$skip);
                $noseg = 0;
            }
        }

        if ($noseg) {
            # No segregating sites in this region, and there is no way to
            # convey information about a homozygous region to msmc.
            next;
        }

        # truncate sequence to match what msmc sees
        splice(@seq, $lastpos+1);

        print $pfh ">$i";
        for my $j (0 .. $#seq) {
            print ($pfh "\n") if ($j % 60 == 0);
            print ($pfh $seq[$j]?'K':'T');
        }
        print $pfh "\n";

        # msmc
        open(my $mfh, ">", "$prefix/msmc.sample$sample.$i")
            or die "cannot open $prefix/msmc.sample$sample.$i for writing: $!";
        my $last = 0;
        for my $x (0 .. $#sites) {
            if (substr($h1,$x,1) != substr($h2,$x,1)) {
                my $pos = $sites[$x] - $from;
                if ($pos == $last) {
                    next;
                }
                my $gap = $pos - $last;
                print $mfh "$i\t$pos\t$gap\t10,01\n";
                $last = $pos;
            }
        }
        close($mfh);
    }

    close($pfh);
}
