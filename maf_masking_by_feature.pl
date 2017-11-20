#!/usr/bin/perl

use strict;

my $usage = "maf_masking_by_feature.pl MAF_file feature_name (FASTA_file GFF_file)xN\n\n This script will take a set of files consisting of a MAF alignment\n and sets of two files: a FASTA file and a GFF annotation file\nof all the sequences that were aligned in said MAF file.\n\n It will then maks the MAF alignment at all positions where the GFF file\nhas a feature_name annotation (typically gene or exon).\n\n";

if ($#ARGV < 3) {
    die $usage;
}

my $MAF = $ARGV[0];
my $feature = $ARGV[1];
my @files;
my %gff;

my %seq;
my %revcomp;

my $fileN = 2;
while ($fileN < $#ARGV) {
    push (@files, $ARGV[$fileN]);
    $gff{$ARGV[$fileN]} = $ARGV[$fileN + 1];
    $fileN+=2;
}

my $name;
foreach my $file (@files) {
    open (READ, $file) || die "Couldn't open sequence file $file\n";
    while (<READ>) {
	chomp;
	if (/>/) {
	    s/>//;
	    s/\s//g;
	    $name = $_;
	}
	else {
	    $seq{$name} .= uc($_); # We make sure it's all uppercase!
	}
    }
    close READ;

    print STDERR "Done with $file\n";

#Now we mask the features we want

#    print STDERR "Starts with: ", length($seq{$name});

    open (READ, $gff{$file}) || die "Couldn't open gff file $gff{$file}\n";
    while (<READ>) {
	if (/^\#/) {}
	else {
	    my @all = split (/\t/, $_);
	    if ($all[2] eq $feature && $seq{$all[0]}) {
		#print $_;
		if ($all[4] < $all[3]) {die;}
		my $start = $all[3];
		my $length = $all[4] - $all[3] + 1;

		# Long-winded way to do this so that it's clear what I'm doing.
		# It can all be compressed into a single substr() operation but
		# this isn't perl golf.

		# Note that we could alternatively just lowercase the string
		# instead of replacing it with multiple "n" characters if
		# it was a different purpose.

		#my $replacement = ("N"x$length);
		my $replacement = ("n"x$length);
		#print substr($seq{$all[0]},($start-1),$length), "\n";
		substr($seq{$all[0]},($start-1),$length,$replacement);
		#print substr($seq{$all[0]},($start-1),$length), "\n";
	    }
	}
    }
    close READ;

#    print STDERR "Ends with: ", length($seq{$name});

    print STDERR "Done with $gff{$file}\n";
}

#print ">maskedchr1\n$seq{$name}\n";
#die;

# Okay, now we have the whole thing masked. Let's go through the MAF file and
# replace the sequence with the now-masked replacement. If the alignment falls
# on a masked region, then we'll have a bunch of Ns!

# Given the way lastz works, the alignment positions on the - strand are
# different (it RCs the entire sequence), so we'll make a RC of the entire
# fucking set of sequences right now to deal with it.

foreach my $n (keys %seq) {
    print STDERR $n, "\n";
    #$revcomp{$n} = $seq{$n};
    $revcomp{$n} = reverse($seq{$n});
    $revcomp{$n} =~ tr/ACTG/TGAC/; # it's all uppercase except for the "n"
}

open (READ, $MAF) || die "Can't open $MAF to read\n";
while (<READ>) {
    if (/^s\s/) {
	my $seq1;
	my $seq2; # These are the two parts of the alignment

	chomp;

	my @all1 = split(/\s+/, $_);
	my $cur = $all1[1];
	my $pos = 0;
	my $ali = $all1[6];
	my $chromosome;

	if ($all1[4] eq "+") {
	    $chromosome = $seq{$cur};
	    if (length($chromosome) < 100) {die "$cur\n";}
	}
	else {
	    $chromosome = $revcomp{$cur};
	}

	# Now comes the replacement;

	for my $n (0..(length($ali)-1)) {
	    my $curpos = $all1[2] + $pos;
	    unless (substr($ali,$n,1) !~ /[a-zA-Z]/) {
		#print "$cur\t$pos\t", length($ali), "\t";
		my $base = substr($chromosome,$curpos,1) || die;
		#print $base; #print "\n";
		substr($ali,$n,1,$base);
		$pos++;
	    }
	}
	$seq1 = $ali;

	$_ = <READ>;
	chomp;
	my @all2 = split(/\s+/, $_);
	$cur = $all2[1];
	$pos = 0;
	$ali = $all2[6];
	$chromosome = "";

	if ($all2[4] eq "+") {
	    $chromosome = $seq{$cur};
	}
	else {
	    $chromosome = $revcomp{$cur};
	}

        # Now comes the replacement;

	for my $n (0..(length($ali)-1)) {
	    my $curpos = $all2[2] + $pos;
	    unless (substr($ali,$n,1) !~ /[a-zA-Z]/) {
		my $base = substr($chromosome,$curpos,1) || die;
		substr($ali,$n,1,$base);
		$pos++;
	    }
	}
	$seq2 = $ali;

	# And now do the reciprocal

	for my $n (0..(length($seq1) - 1)) {
	    if (substr($seq1,$n,1) eq "n" || substr($seq2,$n,1) eq "n") {
		if (substr($seq1,$n,1) ne "-") {
		    substr($seq1,$n,1,"n");
		}
		if (substr($seq2,$n,1) ne "-") {
		    substr($seq2,$n,1,"n");
		}
	    }
	}
	# OK, now print the lines;

	for my $n (0..5) {
	    $all1[$n] =~ s/\s+//g;
	    print $all1[$n], " ";
	}
#	$seq1 =~ s/n/N/g;
	print $seq1, "\n";
	for my $n (0..5) {
	    $all2[$n] =~ s/\s+//g;
	    print $all2[$n], " ";
	}
#	$seq2 =~ s/n/N/g;
	print $seq2, "\n";
    }
    else {
	print;
    }    
}
close READ;
