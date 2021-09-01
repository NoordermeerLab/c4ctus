#!/usr/bin/env perl

use strict;
use warnings;


open(IN, "<" . $ARGV[0]) or die "Impossible to open input file\n";
open(Outfwd, ">" . $ARGV[1]) or die "Impossible to open out fwd file\n";
open(Outrev, ">" . $ARGV[2]) or die "Impossible to open out fwd file\n";

while (my $s = <IN>) {
	chomp $s;
	my @a = split('\t', $s);
	
	if ($a[5] eq "+") {
		print Outfwd join("\t", $a[0], $a[1], $a[2], 1/$a[4]), "\n";
	} elsif ($a[5] eq "-") {
		print Outrev join("\t", $a[0], $a[1], $a[2], 1/$a[4]), "\n";
	}
	# -1 => for 0-based coordinates 
}