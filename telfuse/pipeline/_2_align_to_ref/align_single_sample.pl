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
#or die("Error in command line arguments\n");
or die("Usage: perl $0 --genome GENOME --threads INT <fastq1> <fastq2> <label>\n");

# my $filelist = $ARGV[0];
# my $index = $ARGV[1];

if(scalar(@ARGV) != 3){
	die("Usage: perl $0 --genome GENOME --threads INT <fastq1> <fastq2> <label>\n");
}


my $fastq1 	= $ARGV[0];
my $fastq2 	= $ARGV[1];
my $label  	= $ARGV[2];


# Align using bwa-mem
system("bwa mem -t $threads $genome $fastq1 $fastq2 | samtools view -@ 4 -bS - > $label.bam");
#system("samtools sort -@ 16 -m 4G $label.bam > $label.sorted.bam");
system("samtools sort -@ 4 $label.bam > $label.sorted.bam");
system("samtools index -@ 4 $label.sorted.bam");
system("rm $label.bam");

