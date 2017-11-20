#!/usr/bin/perl

$usage = "merge_and_replace.pl reference alternative MAFfile Identity(0-1)\n";

unless ($#ARGV == 3) {
    die $usage;
}

my %refseq;
my %altseq;

my $name;

open (READ, $ARGV[0]);
while (<READ>) {
    if (/>/) {
	s/\>//;
        my @x = split;
	$name = $x[0];
    }
    else {
	chomp;
	$refseq{$name} .= $_;
    }
}
close READ;

print STDERR "Done reading\n";

open (READ, $ARGV[1]);
while (<READ>) {
    if (s/^>//) {
        my @x = split;
        $name = $x[0];
    }   
    else {
        chomp;
        $altseq{$name} .= $_;
    }   
}

# We no longer uppercase anything! We'll use lowercase to make sure we don't
# recursively replace things into mismatches. Now, with lowercase letters,
# which we use to indicate something that's been replaced already, we keep
# doing replacements only in parts that haven't had them before. This way we
# can go recursively down to ~90 identity and not get shit put in good
# contigs.

my $counter = 0;
my $id = $ARGV[3];

open (READ, $ARGV[2]);
 LOOP: while (<READ>) {
     $counter++;
     if (/^s\s/) {
	 chomp;
	 my @x = split;
	 my $sec = <READ>;
	 chomp($sec);
	 my @y = split(/\s+/, $sec);
	 
	 $counter++;
	 
	 if ($x[3] < 500) {
	     next LOOP;
	 }
	 
	 if ($x[6] eq $y[6]) {
	     next LOOP;
	 }
	 
	 # First, fill gaps on the alternative with the reference.
	 
	 my $total_matches;
	 for my $i (0..length($x[6])) {
	     if (substr($y[6],$i,1) =~ /[nN]/) {
		 $replace = substr($x[6],$i,1);
		 substr($y[6],$i,1,$replace);
	     }
	     if (substr($y[6],$i,1) eq substr($x[6],$i,1)) {
		 $total_matches++;
	     }
	 }

	 
	 if ($total_matches/length($x[6]) < $id) {
	     next LOOP;
	 }
	 # Now, we replace the original sequence with the alternative.
	 # to avoid problems with the regexp length, we'll search for
	 # the first 500 bp, and if it occurrs more than once, we increase
	 # it by 100 bases until we only have one match or we reach 5 kb,
	 # whatever happens first.
	 #
	 # Actually, we don't need to do that. The MAF has the correct
	 # alignment positions, curiously enough.
	 
	 #my $target = -1;
	 #my $initial = 400;
	 #my $stay = 1;
	 
	 #while ($stay) {
	 #    my $matcher = substr($ref{$x[1]},0,$initial);
	 #    my $total = () = ($ref{$x[1]} =~ /$matcher/gi);
	 #    if ($total > 1) {
	 #        $initial += 100;
	 #    }
	 #    elsif ($total == 1) {
	 #        $target = $.
	 #        $stay = 0;
	 
	 # Now, we remove all the gaps from the second sequence:
	 
	 $y[6] =~ s/\-//g;
	 $y[6] = uc($y[6]);
	 
	 # And now we replace with a substr!
	 # First let's verify that the reference and substring are
	 # identical (in the unlikely case we get overlapping matches)
	 $x[6] =~ s/\-//g;
	 $x[6] = uc($x[6]);
	 
	 my $original = substr($refseq{$x[1]}, $x[2], length($x[6]));
	 if ($original eq $x[6]) {
	     # We lowercase the y[6] so that we know where a replacement happened
	     $y[6] = lc($y[6]);
	     substr($refseq{$x[1]}, $x[2], $x[3], $y[6]);
	     print STDERR ".";
	 }
	 else {
	     print STDERR "X";
#	    print substr($x[6], 0, 100), "\t", length(substr($refseq{$x[1]}, $x[2], $x[3])), "\n";
#	    print substr($refseq{$x[1]}, $x[2], 100), "\t", length($x[6]), "\n";;
#	     print $x[6], "\n", $refseq{$x[1]}, "\n"; die;
#	    die;
	 }
     }
}

close READ;

foreach my $n (keys %refseq) {
    print ">$n\n$refseq{$n}\n";
}

print STDERR "\n";









