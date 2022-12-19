#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;


my $inputfile 	= $ARGV[0];
my $label 	= $ARGV[1];
if (scalar @ARGV != 2){
	die "Usage: $0 <input_fasta> <output_label>\n";
}

my $repeat_cutoff = 4;
my $dirname 	= dirname(__FILE__);
my $repeat_count_file = $label . ".repeatcounts";
my $repeat_count_file_filtered =  $label . ".repeatcounts.filtered";
my $repeat_count_file_filtered_readname =  $label . ".repeatcounts.filtered.readname";
my $extracted_fasta = $label . ".telomeric.fasta";

# Check if file is a bam or a fasta file
if($inputfile =~ /\.bam$/){
	# Process bam file
	my $count_repeat_bam = $dirname . "/" . "count_motifs_on_read.bam.pl";
	system("perl $count_repeat_bam $inputfile > $repeat_count_file");
}
else{
	# Process fasta file
	my $count_repeat_fasta = $dirname . "/" . "count_motifs_on_read.fasta.pl";
	system("perl $count_repeat_fasta $inputfile > $repeat_count_file");
}

# Filter out the reads above a cutoff
#system("awk -F \"\\t\" \'{if(\$2>=10) print}\' $repeat_count_file \> $repeat_count_file_filtered");
open(my $REPEATCOUNT, $repeat_count_file) || die $!;
open(my $REPEATCOUNTFILTERED, ">", $repeat_count_file_filtered);
while(my $line = <$REPEATCOUNT>){
	chomp($line);
	my @lineArr = split("\t", $line);
	if($lineArr[1] >= $repeat_cutoff){
		print $REPEATCOUNTFILTERED join("\t", @lineArr) . "\n";
	}
}
close($REPEATCOUNT);
close($REPEATCOUNTFILTERED);

# Get readnames required
system("cut -f1 $repeat_count_file_filtered > $repeat_count_file_filtered_readname");


# Extract telomeric fasta
system("perl $dirname/extract_reads_based_on_readnames.fasta.pl $repeat_count_file_filtered_readname $inputfile > $extracted_fasta");


# Compress fasta file
system("gzip $extracted_fasta");



