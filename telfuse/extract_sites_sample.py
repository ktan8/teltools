import sys
import os
import argparse
from pipeline._1_extract_telomeric.extract_telomeric_reads_from_fastq import extract_telomeric_fastq
from pipeline._1_extract_telomeric.extract_telomeric_reads_from_bam import extract_telomeric_bam
#from 4_identify_softclip_sites.extract_softclipped_sites import extract_softclipped_sites

script_path = os.path.dirname(os.path.realpath(sys.argv[0]))

def main(bamfile, genome, label):
	#####################################
	# 1. Extract telomeric reads for analysis
	#####################################
	# Extract telomeric reads (fastq)
	#extract_telomeric_fastq(fastq1, fastq2, label)

	# Extract telomeric reads (bamfile)
	########extract_telomeric_bam(bamfile, label)

	#####################################
	# 2. Align extracted reads to reference
	#####################################
	fastq1 = label + ".fq1.gz"
	fastq2 = label + ".fq2.gz"
	telomeric_aligned_bam_label = label + ".telomeric.aligned_ref"
	align_sample_script = script_path + "/pipeline/_2_align_to_ref/align_single_sample.pl"
	os.system("perl %s --genome %s %s %s %s" %(align_sample_script, genome, fastq1, fastq2, telomeric_aligned_bam_label))

	#print("perl pipeline/_2_align_to_ref/align_single_sample.pl --genome %s %s %s %s" %(genome, fastq1, fastq2, telomeric_aligned_bam_label))

	#####################################
	# 3. Extract softclipped mappings
	#####################################
	telomeric_aligned_bam = label + ".telomeric.aligned_ref.sorted.bam"
	telomeric_aligned_softclipped_bam = label + ".telomeric.aligned_ref.softclipped.bam"
	extract_softclipped_mapping_script = script_path + "/pipeline/_3_extract_softclip_reads/extract_soft_clipped_mappings.pl"

	os.system("perl %s %s | samtools view -b > %s" %(extract_softclipped_mapping_script, telomeric_aligned_bam, telomeric_aligned_softclipped_bam))

	#####################################
	# 4. Identify sites from softclipped bam
	#####################################
	# python ~/code/FuseTect/identify_intra_telomeric_sites/extract_softclipped_sites.py {} '>' {/.}.analysis.txt
	telomeric_aligned_softclipped_sites = label + ".telomeric.aligned_ref.softclipped.sites"
	extract_softclipped_sites_script = script_path + "/pipeline/_4_identify_softclip_sites/extract_softclipped_sites.py"
	os.system("python %s %s > %s" %(extract_softclipped_sites_script, telomeric_aligned_softclipped_bam, telomeric_aligned_softclipped_sites))


	####################################
	# 5. Aggregate softclipped sites in single sample
	####################################
	telomeric_aligned_softclipped_sites_aggregated = label + ".telomeric.aligned_ref.softclipped.sites.aggregated"
	softclipped_analysis_script = script_path + "/pipeline/_5_aggregate_softclip_reads/softclipped_analysis.py"
	os.system("python %s %s > %s" %(softclipped_analysis_script, telomeric_aligned_softclipped_sites, telomeric_aligned_softclipped_sites_aggregated))








# 1. Extract telomeric short-reads
# 'extract_telomeric_reads_from_fastq.py'
# - Extract telomeric reads
# - uses 2x TTAGGG as a heuristic for search

# 2. Align extracted telomeric reads back to reference genome

# `nohup bash -c "seq 0 195 | parallel -j 15 perl /homes6/kartong/code/helper/alignment/process_alignments.pl --genome /meyersonlab/kartong/genome/Homo_sapiens_assembly38.fasta /meyersonlab-archive/maxgj/tcga_analysis/parallel.files/tcga_hg38_alignment_parallel.txt {}" &`

# 3. Extract softclipped reads

# ls ../bamfile/*.bam | parallel -j 12 perl /homes6/kartong/code/FuseTect/identify_intra_telomeric_sites/extract_soft_clipped_mappings.pl {} '>' {/.}.softclipped.bam


# 4. Identify softclipped sites

# ls ../bam_softclipped/*.bam | parallel -j 12 python ~/code/FuseTect/identify_intra_telomeric_sites/extract_softclipped_sites.py {} '>' {/.}.analysis.txt

# 5. Aggregate softclipped sites in single sample

# ls ../softclipped/*.analysis | parallel -j 24 python /homes6/kartong/code/FuseTect/identify_intra_telomeric_sites/softclipped_analysis.py {} '>' {/.}.aggregated



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Get candidate ectopic telomeric site from a single sample',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)



	parser.add_argument('bamfile', metavar='bamfile', type=str, nargs=1,
		help='bamfile to analyze')
	parser.add_argument('genome', metavar='genome', type=str, nargs=1,
		help='reference genome to map reads onto')
	parser.add_argument('label', metavar='label', type=str, nargs=1,
		help='output label indicating path/prefix for output files')
	# parser.add_argument('--cutoff', metavar='cutoff', type=float, default=0.35,
	# 				help='cutoff of repeat signal')
	# parser.add_argument('--movave', metavar='movave', type=int, default=50,
 #                                        help='Window size for moving average window')
	# parser.add_argument('--movmed', metavar='movmed', type=int, default=501,
 #                                        help='Window size for moving median window')
	# parser.add_argument('--motif', metavar='motif', type=str, default='TTAGGG',
	# 				help='motif to assess for repeat signal (e.g. TTAGGG)')
	# parser.add_argument('--format', metavar='format', type=str, default='pdf',
	# 				help='output file type (e.g. pdf, png)')
	# parser.add_argument('--folder', metavar='folder', type=str, default=image_dir,
 #                                        help='folder for images')
	# parser.add_argument('--noseq', action='store_true',
	# 				help='do not print sequence of read in output')
	# parser.add_argument('--nofig', action='store_false',
	# 				help='do not plot figures for sequences')
	# #parser.add_argument('--fastq', metavar='fastq', type=str, default='./',
 #        #                                help='folder for images')
	# parser.add_argument('--threads', metavar='threads', type=int, default=12,
 #                                        help='number of threads to use')
	# parser.add_argument('--penalty', action='store_true',
	# 				help='calculate telomere length using penalty based method')
	# parser.add_argument('--penaltyval', metavar='movave', type=int, default=0.5,
	# 				help='penalty value to use for penalty method')
	args = parser.parse_args()


	# main(args.fasta[0], pic_format=args.format, AverageTelomereSignalCutoff=args.cutoff,
	# 	short_motif = args.motif, movingAveWindow=args.movave, medianKernelSize=args.movmed,
	# 	folder=args.folder, threads = args.threads, noseq=args.noseq,
	# 	penalty=args.penalty, penaltyval=args.penaltyval, plot_fig=args.nofig)


	main(args.bamfile[0], args.genome[0], args.label[0])
