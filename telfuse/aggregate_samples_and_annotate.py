import sys
import os
import argparse


folder = sys.argv[1]
label = sys.argv[2]
genome = sys.argv[3]


script_path = os.path.dirname(os.path.realpath(sys.argv[0]))

###################################################
# 6. Aggregate multiple samples into a single file
###################################################
aggregated_file = label + ".aggregated_samples.txt"
os.system("python %s/pipeline/_6_aggregate_samples/aggregate_samples.py %s > %s" %(script_path, folder, aggregated_file))


###################################################
# 7. Generate a panel of normal
###################################################
PON_file = label + ".aggregated_samples.PON.txt"
os.system("python %s/pipeline/_7_generate_pon/generate_PON.py %s > %s" %(script_path, aggregated_file, PON_file))


###################################################
# 8. Add edge information
###################################################
aggregated_edge_file = label + ".aggregated_samples.edge.txt"
genome_fai = genome + ".fai"
os.system("python %s/pipeline/_8_add_annotation/annotate_edge.py %s %s > %s" %(script_path, aggregated_file, genome_fai, aggregated_edge_file))


###################################################
# 9. Add PON information
###################################################
aggregated_edge_pon_file = label + ".aggregated_samples.edge.pon.txt"
os.system("python %s/pipeline/_8_add_annotation/annotate_PON.py %s %s > %s" %(script_path, aggregated_edge_file, PON_file, aggregated_edge_pon_file))


###################################################
# 10. Classify telomeric sites
###################################################
aggregated_edge_pon_telo_file = label + ".aggregated_samples.edge.pon.teloclassification.txt"
os.system("python %s/pipeline/_8_add_annotation/classify_telomeric_sites.py %s > %s" %(script_path, aggregated_edge_pon_file, aggregated_edge_pon_telo_file))


###################################################
# 11. Filter for sites meeting criteria
###################################################
aggregated_edge_pon_telo_filtered_file = label + ".aggregated_samples.edge.pon.teloclassification.filtered.txt"
os.system("bash %s/pipeline/_9_filter_sites/filtersites.sh %s > %s" %(script_path, aggregated_edge_pon_telo_file, aggregated_edge_pon_telo_filtered_file))

# os.system("awk -F \"\t\" '{if($15 == 0 && $16==1 && $5>=3 && $6>=0.95){print}}\' aggregated_sites.withMapQ.edge.PON.classification.txt") 
# # > aggregated_sites.withMapQ.edge.PON.classification.filtered.txt

# awk -F "\t" '{if($10>=30){print}}' aggregated_sites.withMapQ.edge.PON.classification.filtered.txt > aggregated_sites.withMapQ.edge.PON.classification.filtered.teloMapQ30.txt
# awk -F "\t" '{if($19!="NA"){print}}' aggregated_sites.withMapQ.edge.PON.classification.filtered.teloMapQ30.txt > aggregated_sites.withMapQ.edge.PON.classification.filtered.teloMapQ30.teloseq.txt

###################################################
# 12. Add genome seq of candidate sites
###################################################
aggregated_edge_pon_telo_filtered_genome_file = label + ".aggregated_samples.edge.pon.teloclassification.filtered.genomeseq.txt"
header = ["Sample", "Chromosome", "Position", "Softclipped_orientation",
"Number_of_supporting_reads", "Average_sequence_identity_of_softclipped_sequences", 
"Average_weighted_sequence_identity_of_softclipped_sequences", "Softclipped_consensus sequence",
"Average_position_of_softclipped_site_on_read",	"Average_mapping_quality_of_reads",
"Min_mapping_quality_of_reads",	"Max_mapping_quality_of_reads",	"Number_of_basepairs_from_left_edge_of_chromosome",
"Number_of_basepairs_from_right_edge_of_chromosome", "Site_is_near_edge_of_chromosome?",
"Number of sites also detected in panel of normal", "Total number of reads corresponding in panel of normal",
"First_12bp_of_softclipped_sequence", "Site_classification", "Genomic_sequence_prior_to_telomeric_event",
"Genomic_sequence_prior_to_telomeric_event (rev_complement)", "Genomic_sequence_after_telomeric_event",
"Genomic_sequence_after_telomeric_event (rev_complement)"]

header_str = "\t".join(header)
os.system("echo \'%s\' > %s" %(header_str, aggregated_edge_pon_telo_filtered_genome_file))
os.system("perl %s/pipeline/_8_add_annotation/extract_ref_sequence.v3.pl %s %s >> %s" %(script_path, aggregated_edge_pon_telo_filtered_file, genome, aggregated_edge_pon_telo_filtered_genome_file))




