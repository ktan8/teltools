#!/usr/bin/perl

use strict;
use warnings;

my $bamfile = $ARGV[0];

open(my $BAMFILE, "-|", "samtools view -h $bamfile") || die $!;

while(my $bamline = <$BAMFILE>){
	if ($bamline =~ /^@/){
		print $bamline;
		next;
	}
	chomp($bamline);
	my @lineArr = split(/\t/, $bamline);
	my $cigar_str = $lineArr[5];
	# Check for a soft clip tag
	if($cigar_str =~ /S/){
		print $bamline . "\n";
	}
}



close($BAMFILE);
