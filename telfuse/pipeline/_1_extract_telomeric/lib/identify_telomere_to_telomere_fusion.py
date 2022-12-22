
from lib.motif_lib import *


def get_telomere_to_telomere_fusion_combinations(motif):
	'''
	Generate the sequence we expect to see
	in the event of a telomere-to-telomere fusion.
	E.g. We should see TTAGGG-CCCTAA, and the full
	array of 36 different combinations.
	'''
	motif_revCom = reverseComplement(motif)
	motif_fwd_combinations = generateMotifCombinations(motif)
	motif_rev_combinations = generateMotifCombinations(motif_revCom)

	# Get all the combinations
	fusion_combinations = []
	for motif_fwd_curr in motif_fwd_combinations:
		for motif_rev_curr in motif_rev_combinations:
			fusion_seq = motif_fwd_curr + motif_rev_curr
			#print fusion_seq
			fusion_combinations.append(fusion_seq)

	return fusion_combinations

def check_telomere_to_telomere_fusion(sequence, motif="TTAGGG"):
	'''
	Identify if it is a telomere-telomere fusion
	based on whether it is a TTAGGG joined to a
	CCCTAA sequence. We look for reads supporting 
	a (TTAGGG)n-(CCCTAA)n junction.
	'''
	fusion_combinations = get_telomere_to_telomere_fusion_combinations(motif)
	for fusion_seq in fusion_combinations:
		if fusion_seq in sequence:
			print fusion_seq
			return 1

	return 0
	


if __name__ == "__main__":
	print get_telomere_to_telomere_fusion_combinations("TTAGGG")
	print check_telomere_to_telomere_fusion("TTAGGGCCCTAA")



