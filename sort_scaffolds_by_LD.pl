#!/usr/bin/perl

my $usage = "sort_scaffolds_by_LD.pl LDfile Chromosomes.fasta unsorted_scaffolds.fasta\n";
if ($#ARGV != 2) {
    die $usage;
}

open (READ, $ARGV[1]);
while (<READ>) {
    chomp;
    if (s/>//) {
	$name = $_;
    }
    else {
	$chr{$name} .= $_;
    }
}
close READ;
print STDERR "Done part 1!\n";

open (READ, $ARGV[2]);
while (<READ>) {
    chomp;
    if (s/>//) {
	$name = $_;
    }
    else {
	$contig{$name} .= $_;
    }
}
close READ;

print STDERR "Done part 2!\n";

open (READ, $ARGV[0]);
my $none = <READ>;
my $line = 1;
while (<READ>) {
    $line++;
    chomp;
    s/\"//g;
    my ($pos, $ori, $new) = split;
    my $chunk = substr($chr{$ori}, ($pos - 100), 200);
    $chunk =~ s/^N+//;
    $chunk =~ s/N+$//;
    my $revchunk = $chunk;
    $revchunk = reverse($revchunk);
    $revchunk =~ tr/ACTG/TGAC/;
    foreach my $n (keys %contig) {
	if ($contig{$n} =~ /$chunk/ || $contig{$n} =~ /$revchunk/) {
	    print ">$n\t$new\t$line\n";
	}
    }
}
close READ;
