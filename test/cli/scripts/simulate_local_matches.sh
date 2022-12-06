#!/usr/bin/env bash
set -e

BINARY_DIR="../../lib/raptor_data_simulation/build/bin"
ERROR_RATE=$1
MATCH_COUNT=$2
MIN_LEN=$3
MAX_LEN=$4

echo "Sampling $MATCH_COUNT local matches between $MIN_LEN and $MAX_LEN bp with an error rate of $ERROR_RATE"

match_dir=local_matches_e$ERROR_RATE
mkdir -p $match_dir
$BINARY_DIR/generate_local_matches \
	--output $match_dir \
	--max-error-rate $ERROR_RATE \
	--num-matches $MATCH_COUNT \
	--min-match-length $MIN_LEN \
	--max-match-length $MAX_LEN \
	--verbose-ids \
	ref.fasta &> /dev/null

mv $match_dir/ref.fastq local_matches/e$ERROR_RATE.fastq
rm -r $match_dir
