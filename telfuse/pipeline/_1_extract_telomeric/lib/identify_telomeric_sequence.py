from .motif_lib import generateMotifCombinations
from .motif_lib import reverseComplement

def check_motif_present(sequence, motif):
	'''
	Check if a sequence has a motif, or
	a rearrangement of it.
	'''
	motif_combinations = generateMotifCombinations(motif)

	for motif in motif_combinations:
		if motif in sequence:
			return 1

	# Return a false if no motif was found
	# across all the combinations
	return 0


def check_motif_present_unstranded(sequence, motif):
	'''
	Check if a motif is present on either the forward
	or the reverse strand.
	'''
	motif_revCom = reverseComplement(motif)
	if check_motif_present(sequence, motif) or check_motif_present(sequence, motif_revCom):
		return 1
	else:
		return 0


def check_telomeric_sequence(sequence, telomere_motif = "TTAGGG", repeat_count=2):
	'''
	Check if a sequence is possibly a telomeric
	sequence. We require it to have two consecutive
	copies of the telomeric repeat sequence. (i.e. 12mer)
	'''

	#long_motif = telomere_motif + telomere_motif
	long_motif = repeat_count * telomere_motif

	if check_motif_present_unstranded(sequence, long_motif):
		return 1
	else:
		return 0



