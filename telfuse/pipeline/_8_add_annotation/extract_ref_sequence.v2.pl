#!/usr/bin/perl

use strict;
use warnings;

my $candidate_list = $ARGV[0];
my $genome 	= $ARGV[1];
my $seq_len 	= 30;

open(my $FILE, $candidate_list) || die $!;

while(my $line = <$FILE>){
	chomp($line);
	my @lineArr = split(/\t/, $line);
	my $chr 	= $lineArr[1];
	my $posn 	= $lineArr[2];
	my $softclipped_direction = $lineArr[3];
	
	my $event_genome_seq = "NA";
	my $event_genome_revcom_seq = "NA";	
	my $lost_genome_seq = "NA";
	my $lost_genome_revcom_seq = "NA";	


	# When softclipped sequence is on the right
	if($softclipped_direction eq "right"){
		my $event_start = $posn - $seq_len;
		my $event_end   = $posn - 1;
		my $event_coord = $chr . ":" . $event_start . "-" . $event_end;

		my $lost_start = $posn;
		my $lost_end   = $posn + $seq_len - 1;
		my $lost_coord = $chr . ":" . $lost_start . "-" . $lost_end;

		my $event_fasta = `samtools faidx -n 1000 -i $genome $event_coord`;
		my ($event_label, $event_sequence) = split(/\n/, $event_fasta);
		$event_genome_seq = $event_sequence;
		$event_genome_revcom_seq = reverse_complement($event_sequence);

		my $lost_fasta = `samtools faidx -n 1000 -i $genome $lost_coord`;
		my ($lost_label, $lost_sequence) = split(/\n/, $lost_fasta);
		$lost_genome_seq = $lost_sequence;
                $lost_genome_revcom_seq = reverse_complement($lost_sequence);
	
	}
	# When softclipped sequence is on the left of the
	# genome, but we will reorient it to be on the right
	# when we are reporting the sequence
	elsif($softclipped_direction eq "left"){
		my $event_start = $posn;
		my $event_end = $posn + $seq_len -1;
		my $event_coord = $chr . ":" . $event_start . "-" . $event_end;
		
		my $lost_start = $posn - $seq_len;
		my $lost_end   = $posn - 1;
		my $lost_coord = $chr . ":" . $lost_start . "-" . $lost_end;

		my $event_fasta = `samtools faidx -n 1000 -i $genome $event_coord`;
		my ($event_label, $event_sequence) = split(/\n/, $event_fasta);
		$event_genome_seq = reverse_complement($event_sequence);
		$event_genome_revcom_seq = $event_sequence;	
		
		my $lost_fasta = `samtools faidx -n 1000 -i $genome $lost_coord`;
		my ($lost_label, $lost_sequence) = split(/\n/, $lost_fasta);
		$lost_genome_seq = reverse_complement($lost_sequence);
		$lost_genome_revcom_seq = $lost_sequence;
	}
		
	push(@lineArr, $event_genome_seq, $event_genome_revcom_seq, $lost_genome_seq, $lost_genome_revcom_seq);
	print join("\t", @lineArr) . "\n";
}


close($FILE);


sub reverse_complement{
	my $dna = shift @_;        

	# Reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
	$revcomp =~ tr/atcgnATCGN/tagcnTAGCN/;
        
	return $revcomp;
}
