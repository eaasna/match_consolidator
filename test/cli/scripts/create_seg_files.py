import pandas as pd
import sys
import os
from Bio import SeqIO

#------------ INPUT ------------ 
#ref_file = "../100kb/ref_rep0.fasta"
ref_file = snakemake.input.ref
#seg_file = "../100kb/split/ref_seg_rep0.txt"
seg_file = snakemake.input.seg_meta
bins = snakemake.config["ibf_bins"]

#------------ OUTPUT ------------
outfile_prefix = snakemake.params.out_prefix
meta_outfile = snakemake.output.meta

# read segment metadata
segments = pd.read_csv(seg_file, header=None, sep='\t')
segments.columns = ["bin_id", "ref_id", "start", "length"]
reference_segments = list(SeqIO.parse(ref_file, "fasta"))

with open(meta_outfile, 'w') as m:
    # assuming a single reference sequence
    for index, row in segments.iterrows():
        bin_id = str(row["bin_id"]).zfill(len(str(bins)))
        ref_id = row["ref_id"]
        start = row["start"]
        length = row["length"]
        fasta_out = outfile_prefix + str(bin_id) + ".fasta"
        m.write(fasta_out + '\n')
        with open(fasta_out, 'w') as f:
            f.write('>' + str(bin_id) + '\t' + "ref=" + str(ref_id) + ",start=" + str(start) + ",length=" + str(length) + '\n')
            f.write(str(reference_segments[0].seq[start:start+length]))
        f.close()
m.close()

