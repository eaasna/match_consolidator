def get_float_er(wildcards):
        if (wildcards.er=="0"):
                # minimum allowed error rate should be 1e-7
                a = 1e-5          # but 1e-7 and 1e-6 lead to invalid pointer error
                return f'{a:.5f}' # supress scientific notation 
        return float(wildcards.er)

rule stellar:
        input:
                ref = "ref.fasta",
                query = "query/with_insertions_e{er}.fasta"
        output: 
                "stellar/e{er}.gff"
        params:
                e = get_float_er
        benchmark:
                "benchmarks/stellar_e{er}.txt"
        shell:
                "stellar --verbose {input.ref} {input.query} -e {params.e} -l {min_len} -a dna -o {output}"

rule dream_stellar:
	input:
		ref = "ref_{bin}.fasta",
		query = "query/with_insertions_e{er}.fasta"
	output:
		"dream_stellar/b{bin}_e{er}.gff"
	params:
		e = get_float_er
	benchmark:
		"benchmarks/dream_stellar_b{bin}_e{er}.txt"
	shell:
		"stellar --verbose {input.ref} {input.query} -e {params.e} -l {min_len} -a dna -o {output}"

rule concat_dream_stellar:
	input:
		expand("dream_stellar/b{bin}_e{{er}}.gff", bin = bin_list)
	output:
		"dream_stellar/joined_e{er}.gff"
	shell:
		"cat {input} > {output}"

# TODO: 
# Merge adjacent alignments when overlapping two segments
# Pick alignments based on search scheme (single best ...) 
rule consolidate_dream_stellar:
	input:
		"dream_stellar/joined_e{er}.gff"
	output:
		"dream_stellar/e{er}.gff"
	shell:
		"consolidate -i {input} -o {output}"

