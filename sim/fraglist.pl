#!/usr/bin/env perl
# Print fragmentation sites, from empirical data.

use strict;
use warnings;
use Getopt::Std;
use List::Util qw(sum min shuffle);

my %opts = (f=>"", l=>0, m=>1e6, s=>31415, x=>0);
getopts("f:l:m:s:x:", \%opts);
my ($faidx, $scale, $lambda) = ($opts{f}, $opts{x}, $opts{l});
my (@fraglens, %masterfraglens);

sub usage {
    print STDERR "Usage:\n";
    print STDERR "       $0 [OPTS..] -f frag.fai master.fai\n";
#    print STDERR "       $0 [OPTS..] -x scale -l lambda master.fai\n";
    print STDERR "\n";
    print STDERR "       -s SEED      random seed [$opts{s}]\n";
    print STDERR "       -m MINLEN    minimum fragment length [$opts{m}]\n";
    exit 1;
}

srand(int($opts{s}));

if ($opts{m} < 1) {
    die "Err: -m $opts{m} is out of range";
}

if ($faidx eq "") {
    if ($scale == 0 || $lambda == 0) {
        usage();
    } else {
        die "Err: not yet implemented";
    }
} else {
    if ($scale != 0 || $lambda != 0) {
        die "Err: -f and -x/-l are mutually exclusive";
    }
    open(my $f, $faidx) or die "Err: cannot open $faidx: $!";
    while (<$f>) {
        chomp;
        my @fields = split(/\s/);
        my $len = int($fields[1]);
        if ($len < int($opts{m})) {
            next
        };
        push(@fraglens, $len);
    }
    close($f);
}

while (<>) {
    chomp;
    my @fields = split(/\s/);
    $masterfraglens{$fields[0]} = int($fields[1]);
}

if (sum(@fraglens) > sum(values(%masterfraglens))) {
    die "Err: total length of fragments too great for master list";
}

@fraglens = shuffle(@fraglens);
my $smallest = min(@fraglens);
# sort before shuffling for reproducibility,
# because key order doesn't respect srand seed
my @mfkeys = shuffle(sort(keys(%masterfraglens)));
my %mfoff;
foreach (@mfkeys) { $mfoff{$_} = 0; }

while (1) {
    my $len = shift(@fraglens);
    my ($k, $ki);
    do {
        $ki = int(rand($#mfkeys));
        $k = $mfkeys[$ki];
        if (not defined $k) {
            die "undefined k: $ki, $#mfkeys";
        }
    } while ($mfoff{$k}+$len >= $masterfraglens{$k});

    my $from = $mfoff{$k}+1;
    my $to = $from+$len;
    $mfoff{$k} = $to;

    print "$k\t$from\t$to\n";

    last if not $#fraglens;

    if ($len == $smallest) {
        $smallest = min(@fraglens);
    }

    if ($mfoff{$k}+$smallest >= $masterfraglens{$k}) {
        splice @mfkeys, $ki, 1;
        delete $mfoff{$k};
    }
}
