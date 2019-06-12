#!/usr/bin/perl

my @total;

# @total is three vaules: Self-self identicals, self-self too far away,
# self-self in tandem.

open (READ, $ARGV[0]);
LOOP: while (<READ>) {
    if (/^a\s/) {
	my $line_ref = <READ>;
	my @ref = split(/\s+/, $line_ref);
	my $line_query = <READ>;
	my @query = split (/\s+/, $line_query);

	# Check if it's a self alignment

	my $selfie = 0;
	if ($ref[1] eq $query[1]) {
	    $selfie = 1;
	}
	unless ($selfie) {
	    next LOOP;
	}

	# Check if they're in tandem, 25 bp gap HARCODED
	# change here if you want to:
	my $gapsize = 25;

	if ($ref[2] == $query[2]) {
	    $total[0]++;
	    next LOOP;
	}
	
	my $tandem_gap = $query[2] - ($ref[2] + length($ref[6]));
	
	if (abs($tandem_gap) <= $gapsize) {
	    $total[2]++;
	    next LOOP;
	}
	if ($tandem_gap > $gapsize) {
	    $total[1]++;
	    next LOOP;
	}
    }
}

print "
Total self-matches:   $total[0]
Total far-matches:    $total[1]
Total tandem-matches: $total[2]\n";
