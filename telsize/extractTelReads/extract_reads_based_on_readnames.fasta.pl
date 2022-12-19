#!/usr/bin/perl

use strict;
use warnings;

my $readnamefile = $ARGV[0];
my $fastafile = $ARGV[1];


open(my $READNAMES, $readnamefile) || die $!;
my %readname_hash;
while(my $line = <$READNAMES>){
	chomp($line);
	$readname_hash{$line} = 1;
}
close($READNAMES);


open(my $FASTA, "-|", "zcat -f $fastafile") || die $!;
while(my $readname = <$FASTA>){
	chomp($readname);
	my $readname_clean = substr($readname, 1);
	my $sequence = <$FASTA>;
	chomp($sequence);

        # Print reads if name matches
        if(exists($readname_hash{$readname_clean})){
        	print ">" . $readname_clean . "\n";
        	print $sequence . "\n";
        }
}
close($FASTA);

