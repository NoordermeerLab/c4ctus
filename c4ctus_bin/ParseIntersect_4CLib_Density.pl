#!/usr/bin/env perl



use strict;
use warnings;




my $FileIN = $ARGV[0];
open(IN, "<" . $FileIN) or die "Impossible to open FileIN \n";



my @AllScore = ();

#Store First Line
my $line = <IN>;
chomp $line;
my @A = split("\t",$line);
my $StoredChr    = $A[0];
my $StoredStart  = $A[1];
my $StoredEnd    = $A[2];
my $StoredInfo   = $A[3];
my $StoredScore  = $A[4];
my $Score        = $A[4];


# Process the rest of the file
while ($line = <IN>) {
	chomp $line;
	@A = split("\t",$line);
	#Current info
	my $CurrentChr    = $A[0];
	my $CurrentStart  = $A[1];
	my $CurrentEnd    = $A[2];
	my $CurrentInfo   = $A[3];
	my $CurrentScore  = $A[4];
	push (@AllScore, $CurrentScore);
	
	
	if ($CurrentChr eq $StoredChr && $CurrentStart == $StoredStart && $CurrentEnd == $StoredEnd && $CurrentInfo eq $StoredInfo) {
		#Position and Fragment type are the same than before => duplicated entry
		$Score = $Score + $CurrentScore;
	}
	else
	{
		print $StoredChr,"\t",$StoredStart,"\t",$StoredEnd,"\t",$StoredInfo,"\t",$Score/($StoredEnd-$StoredStart),"\n";
		$StoredChr    = $CurrentChr;
		$StoredStart  = $CurrentStart;
		$StoredEnd    = $CurrentEnd;
		$StoredInfo   = $CurrentInfo;
		$StoredScore  = $CurrentScore;
		$Score        = $CurrentScore;
		@AllScore = ();
	}
}
print $StoredChr,"\t",$StoredStart,"\t",$StoredEnd,"\t",$StoredInfo,"\t",$Score/($StoredEnd-$StoredStart),"\n";









# 
# #Store First Line
# my $line = <IN>;
# chomp $line;
# my @A = split("\t",$line);
# my $PreviousChr    = $A[0];
# my $PreviousStart  = $A[1];
# my $PreviousEnd    = $A[2];
# my $PreviousInfo   = $A[3];
# #$PreviousInfo =~ /^type=((start|end)Segment)/;
# #my $PreviousType   = $1;
# my $PreviousScore  = $A[4];
# my $Score          = $A[4];
# 
# 
# #Process the rest of the file
# while ($line = <IN>) {
# 	#Store current position, fragment type and score
# 	chomp $line;
# 	@A = split("\t",$line);
# 	my $CurrentChr    = $A[0];
# 	my $CurrentStart  = $A[1];
# 	my $CurrentEnd    = $A[2];
# 	my $CurrentInfo   = $A[3];
# 	#$CurrentInfo =~ /^type=((start|end)Segment)/;
# 	#my $CurrentType   = $1;
# 	my $CurrentScore  = $A[4];
# 	
# 	
# 	if ($CurrentChr eq $PreviousChr && $CurrentStart == $PreviousStart && $CurrentEnd == $PreviousEnd && $CurrentInfo eq $PreviousInfo) {
# 		#Position and Fragment type are the same than before
# 		$Score = $Score + $CurrentScore;
# 	} 
# 	else 
# 	{
# 		#Position and Fragment type are the different than before
# 		if ($Score != 0) {
# 			print $PreviousChr,"\t",$PreviousStart,"\t",$PreviousEnd,"\t",$PreviousInfo,"\t",$PreviousScore,"\n";
# 			print $CurrentChr,"\t",$CurrentStart,"\t",$CurrentEnd,"\t",$CurrentInfo,"\t",$CurrentScore,"\n";
# 			print "---------\n";
# 			print $PreviousChr,"\t",$PreviousStart,"\t",$PreviousEnd,"\t",$PreviousInfo,"\t",$Score/($PreviousEnd-$PreviousStart),"\n\n\n";
# 		}
# 	}
# 	
# 	#store current position for new line
# 	$PreviousChr    = $CurrentChr;
# 	$PreviousStart  = $CurrentStart;
# 	$PreviousEnd    = $CurrentEnd;
# 	$PreviousInfo   = $CurrentInfo;
# 	#$PreviousType   = $CurrentType;
# 	$Score          = $CurrentScore;
# 	$PreviousScore  = $CurrentScore;
# }
# 
# print $PreviousChr,"\t",$PreviousStart,"\t",$PreviousEnd,"\t",$PreviousInfo,"\t",$Score/($PreviousEnd-$PreviousStart),"\n";











