#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;


# Default arguments
my $genome = "~/genome/Homo_sapiens_assembly38.fasta";
my $threads = 40;

# Get arguments
GetOptions ("genome=s"   => \$genome,
	    "threads=i"  => \$threads)
or die("Error in command line arguments\n");

my $filelist = $ARGV[0];
my $index = $ARGV[1];

open(my $INPUT, $filelist) || die $!;


my $currIndex = 0;
while(my $line = <$INPUT>){
	if($currIndex == $index){ 
		chomp($line);
		my @lineArr = split(/\t/, $line);
		my $fastq1 = $lineArr[1];
		my $fastq2 = $lineArr[2];
		my $label = $lineArr[0];
		#system("python /homes6/kartong/code/FuseTect/extract_telomeric_reads_from_fastq.py $fastq1 $fastq2 $label");
		# Align using bwa
		system("bwa mem -t $threads $genome $fastq1 $fastq2 | samtools view -@ 4 -bS - > $label.bam");
		system("samtools sort -@ 16 -m 4G $label.bam > $label.sorted.bam");
		system("samtools index -@ 4 $label.sorted.bam");
		system("rm $label.bam");
	}
	$currIndex += 1;
}

close($INPUT);
