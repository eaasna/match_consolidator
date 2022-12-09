rule compare_output:
	input:
		stellar = "stellar/e{er}.gff",
		dream = "dream_stellar/e{er}.gff"
	output:
		stellar = "evaluate/stellar_e{er}_diff.gff",
		dream = "evaluate/dream_stellar_e{er}_diff.gff"		
	shell:
		"""
		grep -v -f {input.stellar} {input.dream} > {output.dream}
		grep -v -f {input.dream} {input.stellar} > {output.stellar}
		"""

