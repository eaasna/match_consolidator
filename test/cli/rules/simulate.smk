import random
# simulation needs unique seeds otherwise the same sequence is simulated
def get_seed(wildcards):
        return random.randint(0, 1e6)

rule simulate_sequences:
	output:
		ref = "ref.fasta",
		prelim = temp("query/query.fasta"),
		query = temp("query/one_line.fasta")
	params:
		#ref_seed = get_seed,
		#query_seed = get_seed
		ref_seed = 952267,
		query_seed = 541683
	shell:      
		"""
		echo {params.ref_seed}
		echo {params.query_seed}
		../cli/scripts/simulate_sequences.sh {ref_len} {query_len} {params.ref_seed} {params.query_seed}
		"""

rule split_reference:
	input:
		"ref.fasta"
	output:
		bins = expand("ref_{bin}.fasta", bin = bin_list),
		ref_meta = temp("reference_meta.txt"),
		seg_meta = temp("segment_meta.txt"),
		seg_paths = temp("seg_files.txt")
	shell:
		"valik split {input} --reference-output {output.ref_meta} --segment-output {output.seg_meta} --overlap {seg_overlap} --bins {bin_nr}"

rule simulate_matches:
	input:
		ref = "ref.fasta"
	output:
		matches = temp("local_matches/e{er}.fastq")
	params:
		#seed = get_seed
		seed = 391244
	shell:      
		"""
		echo {params.seed}
		../cli/scripts/simulate_local_matches.sh {wildcards.er} {matches} {min_len} {max_len} {params.seed}
		"""

rule insert_matches:
	input:
		query = "query/one_line.fasta",
		matches = "local_matches/e{er}.fastq"
	output:     
		query = "query/with_insertions_e{er}.fasta",
		ground_truth = temp("ground_truth/e{er}.tsv")
	params:     
		#seed = get_seed
		seed = 92436
	script:     
		"../scripts/insert_local_matches.py"

