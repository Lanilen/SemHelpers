#!/usr/bin/perl
#use strict;
my $pos = 0; # Position on genome
my $bases = 50; # Required number of bases with good matching before an insert
my $buffer = 0; # Number of good bases we've found in a row
my $postbuffer = 0; # Number of good bases after a gap
my $threshold = 3/4; # % of bases matching the ref.
my $reference = "";
my $max_mismatches = 1; # This is the number of non-match features that we will
                        # allow when checking the forward/backwards good match
                        # window around discrepancies.

my $state = 0; # States: 
               # 0 = reading matches.
               # 1 = reading inserts.
               # 2 = verifying $bases bases of matches after 1.

my @gap_data = qw(); # Here we'll temporarily store the data that will fill
                     # the gaps.

my $ref_name;

open (READ, $ARGV[0]);
while (<READ>) {
    chomp;
    if (/>/) {
	$ref_name = $_;
    }
    else {
	#chomp;
	$reference .= $_;
    }
}
close READ;

my @alignment; # We're going to get the whole alignment in an array for added
               # convenience. We have to do a lot of forward/backwards checking
               # So it is easier that way.

open (READ, $ARGV[1]);
while (<READ>) {
    push (@alignment, $_);
}
close READ;

BIGLOOP: for my $n (0..$#alignment) {
    $alignment[$n] =~ s/\^~//g;
    my ($name, $new_pos, $ref_base, $depth, $matches,$qual) = split (/\s/, $alignment[$n]);
    $matches = uc($matches);
    $matches =~ s/\$//g;

    # We are going to remove the -XXblahblahblah part, because we don't need
    # it. If a read is missing something in the ref, it'll show up as an 
    # asterisk! So there!

    # Apparently we do! As the rules were laid at the time of this edit, if
    # the _end_ of a sequence is misaligned, the rules weren't checking, and
    # whole chunks of good assembly will be removed for no good reason.

    # So don't use this method, deletions have to be handled the same way as
    # insertions.

#    while ($matches =~ /.\-(\d+)/) {
#	$matches = ${^PREMATCH};
#	$matches .= substr(${^MATCH},0,1);
#	$matches .= substr(${^POSTMATCH},$1);
#    }

    if ($new_pos - $pos == 1) {
	$pos = $new_pos;
	# The position in the alignment will skip parts with 0 coverage, so
	# we do have to check that two consecutive lines are off by only 1 base
	if (length($matches) == $depth) {
	    # We haven't found any inserts in the matching bases, all we have
	    # is [.,ACTGactgNn]
	    my $tmp_bases = $matches;
	    my $num_matches = ($tmp_bases =~ s/[\.\,]//g);
	    my $ratio = $num_matches/$depth;
	    if ($ratio > $threshold) {
		# At least $treshold% of reads have a perfect match to ref
		# Do nothing, this is business as usual
		$buffer++;
	    }
	    elsif ((1-$ratio) > $threshold) {
		# This is the opposite. The ratio of matches vs. non-matches
		# is under the threshold (the majority of reads do NOT match
		# the reference).

		# First thing first, if it's an N, we will rewind until there
		# are no Ns. And if there's perfect matching _there_, we go to
		# the end of the Ns, and check if there's an insert or deletion.
		# If there is, jump ahead of deletion, or check right after
		# insertion. If they all match, you're set.

		# The reason for this is that, under typical moleculo aligns
		# the only time you get disagreements are:
		# 1) Gap in reference that is filled.
		# 2) End-of-read shennanigans.
		# 3) miscalled base in ref.

		# 1 is the typical scenario where the reads will transition
		# into a bunch of Ns. We'll agree to replace the Ns with bases
		# If we match well before the Ns, and if the stretch of Ns
		# ends in $bases agreement to reference _After_ the insert/del
		my $num_gaps = () = $matches =~ /\*/gi;
		if (($num_gaps/$depth) > $threshold) {
		    next BIGLOOP;
		}
		my $proceed = &mismatch ($n, $bases, $threshold, \@alignment);

		unless ($proceed > 0) {
		    next BIGLOOP;
		}

		$matches = uc($matches);
		my %base_hash;
		for my $i (split //, $matches) {
		    $base_hash{$i}++;
		}
		my @keys = sort {$base_hash{$b} <=> $base_hash{$a}} keys (%base_hash);
		if ($ref_base =~ /N/) {
		    if ($base_hash{$keys[0]}/$depth > 0) {
			if ($keys[0] =~ /[actgACTG]/) {
			    substr($reference, ($new_pos - 1), 1, $keys[0]);
			    $buffer++;
			}
			elsif ($keys[0] eq "*") {
			    substr($reference, ($new_pos - 1), 1, "-");
			    $buffer++;
			}
		    }
		    else {
			$buffer = 0;
		    }
		}
		elsif ($ref_base =~ /[ACTGactg]/) {
		    # Seems like we have a mismatch, if the ratio is good enough
		    # and we've passed the "proceed" test, we should replace
		    # it!
		    substr($reference, ($new_pos - 1), 1, $keys[0]);
		    $buffer = 0;
		}
		else {
		    $buffer = 0;
		}
	    }
	    else {$buffer = 0;}
	}
	else {
	    # There's an insert here! So we have to compare the different
	    # inserts, and see which one we choose.

	    # There could also be deletions! Let's check which one is most
	    # common.

	    # First, we check that there are only identical matches
	    # (i.e., [.,]+XX, where XX is the number of bases)

	    my $proceed = &mismatch ($n, $bases, $threshold, \@alignment);
	    # first, we check that there is agreement up and downstream from
	    # this spot.

	    unless ($proceed > 0) {
		next BIGLOOP;
	    }
	    my $num_posi = () = $matches =~ /\+/gi;
	    my $num_nega = () = $matches =~ /\-/gi;

	    my $operation = "Insert";
	    if ($num_posi < $num_nega) {
		$operation = "Delete";
	    }

	    my $total_inserts = 0;
	    my @all_inserts;
	    my $base_corrected = 0;
	    while ($matches =~ m/.[\+\-](\d+)/g) {
		$total_inserts++;
		my $base = substr(${^MATCH},0,1);
		my $insertion = substr(${^POSTMATCH},0,$1);
		if ($base =~ /[ACTGactg]/) {
		    if ($ref_base eq "N") {
			if ($operation eq "Insert") {
			    $operation = "InsertD";
			}
			elsif ($operation eq "Delete") {
			    $operation = "DeleteD";
			}
		    }
		    if ($base_corrected == 0) {
			$insertion = "$base$insertion";
			$base_corrected = 1;
		    }
		}
		push (@all_inserts, $insertion);
	    }
	    my $verify = $total_inserts/$depth;
	    if ($verify > $threshold) {
		#print STDERR "$total_inserts\n$depth\n$threshold\n";
		my $long_pos = 0;
		my $longest = $all_inserts[$long_pos];
		if ($#all_inserts > 0) {
		    my $state = 1;
		    for my $n (0..$#all_inserts) {
			$longest = &LCS($all_inserts[$n], $all_inserts[$n - 1]);
		    }
		}
		#print STDERR $longest, "\n";
		my $correction = "$new_pos\t$operation\t$longest";
		print STDERR "$alignment[$n]$correction\n";
		push (@gap_data, $correction);
	    }
	}
    }
    else {
	$pos = $new_pos;
	$state = 0;
	$buffer = 0;
	$postbuffer = 0;
    }
}

for my $n (0..$#gap_data) {
    my $correction = pop @gap_data;
    my @instructions = split (/\t/, $correction);

    if ($instructions[1] eq "Insert") {
	substr($reference,($instructions[0]),0,$instructions[2]);
    }
    elsif ($instructions[1] eq "InsertD") {
	substr($reference,($instructions[0] - 1),1,$instructions[2]);
    }
    elsif ($instructions[1] eq "Delete") {
	substr($reference,($instructions[0]),length($instructions[2]),"");
    }
    elsif ($instructions[1] eq "DeleteD") {
	substr($reference,($instructions[0] - 1),(length($instructions[2])),substr($instructions[2],0,1));
    }
}

$reference =~ s/\-//g;

print $ref_name;
print "_corrected\n";
while ($reference =~ m/.{1,80}/g) {
    print ${^MATCH}, "\n";
}

sub mismatch {
    # This subroutine checks that the reference sequence should be replaced
    # with the mismatch the assembled reads suggest.

    my ($j, $b, $thr, $ali) = @_;
    my $checked = 0;
    my $mismatches = 0;
    my $original_position = $j;

    # Start by checking backwards

    my @original = split (/\s/, $$ali[$j]); # Quick array to store original data

    #my $base_tracker; # Check that we move one base at a time when we jump back
    #$base_tracker = $original[1];

    # Just realized that ^^ is pointless.

    CHECK: while ($checked < $b) {
	$j--;
	my $current_line = $$ali[$j];

	# First, check that we haven't run into an insertion or deletion.
	# If we have, we'll have to skip them

	my ($curr_name, $curr_new_pos, $curr_ref_base, $curr_depth, $curr_matches,$curr_qual) = split (/\s/, $current_line);
	$curr_matches = uc($curr_matches);

	my $gaps = () = $curr_matches =~ /\*/gi;
	if (($gaps/$curr_depth) > $thr) {
	    next CHECK;
	    # This will speed up flying over long gaps
	}

	my $positives = () = $curr_matches =~ /\+/gi;
	my $negatives = () = $curr_matches =~ /\-/gi;

	if (($positives/$curr_depth) > $thr) {
	    my $tmp_bases;
	    while (m/.\+/g) {
		$tmp_bases .= substr(${^MATCH},0,1);
	    }
	    $curr_matches = $tmp_bases;
	    # Done like this, we jump over inserts because we _know_ they are
	    # all over, and what matters is if the plain bases match nicely
	    # with the reference.
	}
	elsif (($negatives/$curr_depth) > $thr) {
	    my $tmp_bases;
	    while (m/.\-/g) {
		$tmp_bases .= substr(${^MATCH},0,1);
	    }
	    $curr_matches = $tmp_bases;
	    # After a myriad gaps we will get to the actual -XXblahblah part
	    # Which we will treat as the +XXblahblah, what matters is the
	    # matching on both ends of gaps and inserts. Badly matched
	    # read ends with long gaps will have ONE end that matches badly, so
	    # bypasing te whole gap is fine.
	}
	if ($curr_ref_base !~ /N/ && ($gaps/$curr_depth) < $thr) {
	    # We are neither matching against an N (which, of course, will
	    # lead to a mismatch we WANT to replace) nor against a gap (which
	    # we fly over, because the important part of gaps is the ends, not
	    # the middle).
	    $checked++;
	    my $curr_num_matches = () = $curr_matches =~ /[\.\,]/gi;
	    my $curr_ratio = $curr_num_matches/$curr_depth;
	    if ($curr_ratio < $thr) {
		$mismatches++;
	    }
	}
    }
    if ($mismatches > 1) {
	#print STDERR $$ali[$original_position], "0_1\n";
	return 0;
    } # One of the ends is already bad! So, don't do the replacement.

    # ALTER NUMBER OF ALLOWED MISMATCHES HERE

    $checked = 0;
    $j = $original_position;
    $mismatches = 0;

    CHECK: while ($checked < $b) {
	$j++; # Now we check the other way! We still do things the same, try
	      # and run along, looking for the number of mismatches
	my $current_line = $$ali[$j];

	# First, check that we haven't run into an insertion or deletion.
	# If we have, we'll have to skip them

	my ($curr_name, $curr_new_pos, $curr_ref_base, $curr_depth, $curr_matches,$curr_qual) = split (/\s/, $current_line);
	$curr_matches = uc($curr_matches);

	my $gaps = () = $curr_matches =~ /\*/gi;
	if ($curr_depth > 0) {
	    if (($gaps/$curr_depth) > $thr) {
		next CHECK;
		# This will speed up flying over long gaps
	    }
	}

	my $positives = () = $curr_matches =~ /\+/gi;
	my $negatives = () = $curr_matches =~ /\-/gi;

	if (($positives/$curr_depth) > $thr) {
	    my $tmp_bases;
	    while (m/.\+/g) {
		$tmp_bases .= substr(${^MATCH},0,1);
	    }
	    $curr_matches = $tmp_bases;
	    # Done like this, we jump over inserts because we _know_ they are
	    # all over, and what matters is if the plain bases match nicely
	    # with the reference.
	}
	elsif (($negatives/$curr_depth) > $thr) {
	    my $tmp_bases;
	    while (m/.\-/g) {
		$tmp_bases .= substr(${^MATCH},0,1);
	    }
	    $curr_matches = $tmp_bases;
	    # When we find a -XXblahblah, it'll be followed as gaps, which we
	    # treat as usual up there. So, just blaze by them!
	}
	if ($curr_ref_base !~ /N/ && ($gaps/$curr_depth) < $thr) {
	    # We are neither matching against an N (which, of course, will
	    # lead to a mismatch we WANT to replace) nor against a gap (which
	    # we fly over, because the important part of gaps is the ends, not
	    # the middle).
	    $checked++;
	    my $curr_num_matches = () = $curr_matches =~ /[\.\,]/gi;
	    my $curr_ratio = $curr_num_matches/$curr_depth;
	    if ($curr_ratio < $thr) {
		$mismatches++;
	    }
	}
    }
    if ($mismatches > 1) {
	#print STDERR $$ali[$original_position], "0_2\n";
	return 0;
    } # One of the ends is already bad! So, don't do the replacement.

    # CHANGE NUMBER OF ALLOWED MISMATCHES HERE!

    return (1);
}
    
	    
	

sub LCS {
    my ($a, $b) = @_;
    my $S = [];   # An array of scores
    my $R = [];   # An array of backtracking arrows
    my $n = length($a);
    my $m = length($b);
	
    # We need to work in letters, not in strings.  This is a simple way
    # to turn a string of letters into an array of letters.
    my @a = split // => $a;
    my @b = split // => $b;

    # These are "constants" which indicate a direction in the backtracking array.
    my $NEITHER     = 0;
    my $UP          = 1;
    my $LEFT        = 2;
    my $UP_AND_LEFT = 3;
    
    # It is important to use <=, not <.  The next two for-loops are initialization
    for(my $ii = 0; $ii <= $n; ++$ii) { 
	$S->[$ii][0] = 0;
	$R->[$ii][0] = $UP;
    }
    
    for(my $jj = 0; $jj <= $m; ++$jj) { 
	$S->[0][$jj] = 0;
	$R->[0][$jj] = $LEFT;
    }
    
    # This is the main dynamic programming loop that computes the score and
    # backtracking arrays.
    for(my $ii = 1; $ii <= $n; ++$ii) {
	for(my $jj = 1; $jj <= $m; ++$jj) { 
	    
	    ($S->[$ii][$jj], $R->[$ii][$jj]) = ($a[$ii-1] eq $b[$jj-1]) ?
		($S->[$ii-1][$jj-1] + 1, $UP_AND_LEFT) :
		($S->[$ii-1][$jj-1] + 0, $NEITHER);
	    
	    ($S->[$ii][$jj], $R->[$ii][$jj]) = $S->[$ii-1][$jj] >= $S->[$ii][$jj] ?
		($S->[$ii-1][$jj], $UP) : 
		($S->[$ii][$jj], $R->[$ii][$jj]);
	    
	    ($S->[$ii][$jj], $R->[$ii][$jj]) = $S->[$ii][$jj-1] >= $S->[$ii][$jj] ?
		($S->[$ii][$jj-1], $LEFT) : 
		($S->[$ii][$jj], $R->[$ii][$jj]);
	}
    }
    
    # The length of the longest substring is $S->[$n][$m]
    $ii = $n; 
    $jj = $m;
    my $lcs = '';
    
    # Trace the backtracking matrix.
    while( $ii > 0 || $jj > 0 ) {
	if( $R->[$ii][$jj] == $UP_AND_LEFT ) {
	    $ii--;
	    $jj--;
	    $lcs = $a[$ii].$lcs;
	}
	
	elsif( $R->[$ii][$jj] == $UP ) {
	    $ii--;
	}
	
	elsif( $R->[$ii][$jj] == $LEFT ) {
	    $jj--;
	}
	
	else {
	    die("Uninitialized arrow at ($ii, $jj): $S->[$ii][$jj] / $R->[$ii][$jj]\n");
	}
    }
    
    return $lcs;
}
