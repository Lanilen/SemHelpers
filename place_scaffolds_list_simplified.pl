#!/usr/bin/perl

open (READ, $ARGV[0]);
while (<READ>) {
    chomp;
    s/>//g;
    @data = split;
    push (@original, $data[0]);
    $ref{$data[0]} = $#data;
    $ref_chromosome{$data[0]} = $data[1];
}
close READ;

print STDERR "loaded Ref\n";

my %refhash;

open (READ, $ARGV[1]);
while (<READ>) {
    chomp;
    s/>//g;
    my @data = split;
    if ($ref{$data[0]} > 0 && length($refhash{$data[1]}{$data[8]}) <= 0) {
	$refhash{$data[1]}{$data[8]} = $data[0];
    }
}
close READ;

print STDERR "placed Ref\n";

my %blasthash;
my %placed;
open (READ, $ARGV[1]);
while (<READ>) {
    chomp;
    s/>//g;
    my @data = split;
    unless ($placed{$data[0]}) {
	my %chr_points = %{$refhash{$data[1]}};
	my $closest = g($data[8], keys %chr_points);
	$closest = $refhash{$data[1]}{$closest};
	$blasthash{$closest}{$data[1]}{$data[8]} = $data[0];
	$placed{$data[0]} = 1;
    }
}
close READ;

print STDERR "placed All\n";

# We have all the hits in a hash! Time to show them the ropes!

my $current_chr = "chr0";

foreach my $n (@original) {
    if ($ref_chromosome{$n} ne $current_chr) {
	$current_chr = $ref_chromosome{$n};
	print ">$current_chr\n";
    }
    my @sorter;
    my @name_sorter;
    my $first_lvl = $blasthash{$n};

    for my $key2 (keys %{ $first_lvl }) {
	my $second_lvl = $first_lvl->{$key2};

	for my $key3 (keys %{ $second_lvl }) {
	    push (@sorter, $key3);
	    push (@name_sorter, $second_lvl->{$key3});
	}
    }

    my @idx = sort { $sorter[$a] <=> $sorter[$b] } 0..$#sorter;
    @sorter = @sorter[@idx];
    @name_sorter = @name_sorter[@idx];

    for my $i (0..$#sorter) {
	print $name_sorter[$i], "\t", $sorter[$i];
	if ($ref{$name_sorter[$i]} > 0) {
	    print "\tMARKER\n";
	}
	else {
	    print "\n";
	}
    }
}



sub g{(sort{abs$a-$_[0]<=>abs$b-$_[0]}@_)[1]}
