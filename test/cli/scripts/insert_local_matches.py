import sys
import os
from Bio import SeqIO
import pandas as pd
import random

# ------- INPUT -------
# 1MB of sequence where local matches will be inserted into
query_file = snakemake.input.query
# Local matches sampled with certain error rate
match_file = snakemake.input.matches

# ------- OUTPUT ------- 
query_with_insertions_file = snakemake.output.query
ground_truth_file = snakemake.output.ground_truth

def check_file_exists(fn):
    if not os.path.exists(fn):
        raise SystemError("Error: File does not exist\n")
        
check_file_exists(match_file)
local_matches = list(SeqIO.parse(match_file, "fastq"))

check_file_exists(query_file)
query = list(SeqIO.parse(query_file, "fasta"))
assert (len(query) == 1),"Query contains more than one sequence, pick which one to insert to."

# Sort random positions
query_len = len(query[0].seq)
random.seed(snakemake.params.seed)
ran_ind_list = random.sample(range(query_len), len(local_matches))
ran_ind_list.sort()
# Shuffle local matches so that they would appear in the query in a random order
random.shuffle(local_matches) 

# Assign random location to each local match
for i in range(len(local_matches)):
    
    local_matches[i].random_location = ran_ind_list[i]

original_query = query[0].seq

id_list = []
position_list = []
length_list = []
reference_position_list = []
reference_length_list = []
errors_list = []

# Insert local matches into query at random positions
query_with_insertions = original_query[0:local_matches[0].random_location]
insertion_length = 0
for i in range(0, len(local_matches)):
    match = local_matches[i]
    
    # Gather ground truth
    id_list.append(match.name)
    
    # converts "l50-3 start_position=935090,length=50,errors=0,reference_id='1',reference_file='ref_0.fasta'"
    # into {'start_position': '935090', 'length': '50', 'errors': '0', reference_id="'reference_id'", reference_file="'ref_0.fasta'"}
    meta_info = dict((key.strip(), value.strip()) for key, value in (element.split('=') for element in match.description.split(" ", 2)[1].split(',')))
    reference_position_list.append(int(meta_info['start_position']))
    reference_length_list.append(int(meta_info['length']))
    errors_list.append(int(meta_info['errors']))

    insertion_position = match.random_location + insertion_length
    position_list.append(insertion_position)    # position in the query
    length_list.append(len(match.seq))
    
    prefix = query_with_insertions[0:insertion_position]
    insertion = match.seq
    
    if (i < (len(local_matches) - 1)):
        postfix = original_query[match.random_location:local_matches[i+1].random_location]
    else:
        # Edge case: insert last local match
        postfix = original_query[local_matches[-1].random_location:]
    
    insertion_length += len(match.seq)
    query_with_insertions = prefix + insertion + postfix 

assert (len(query_with_insertions) - len(original_query)) == insertion_length,"Wrong length of sequence inserted."

# Write insertion query file
with open(query_with_insertions_file, 'w') as f:
    f.write('>1\n')
    f.write(str(query_with_insertions))
    
f.close()

# Create ground truth file
ground_truth = {'id':id_list,
                'position':position_list,
                'length':length_list,
                'reference_position': reference_position_list,
                'reference_length': reference_length_list,
                'errors': errors_list}
ground_truth_df = pd.DataFrame(ground_truth)
ground_truth_df.to_csv(ground_truth_file, index=False, sep='\t')
ground_truth_df.head()
