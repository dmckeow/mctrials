#!/bin/bash
eval "$(conda shell.bash hook)"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
echo "Using ${THREADS} cores"

# Defaults
output_dir="1_MI_MCH"
lib_file="../curated_sequences_NR.fa"
log_dir=".."

# Parse command-line options
while getopts ":g:o:l:L:h" option; do
    case "${option}" in
        g) genome=${OPTARG} ;;
        o) output_dir=${OPTARG} ;;
        l) lib_file=${OPTARG} ;;
        L) log_dir=${OPTARG} ;;
        h) echo "Usage: bash script.sh -g <genome.fa> [-o <output_dir>] [-l <library.fa>] [-L <output directory for logfile>]"; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG"; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument."; exit 1 ;;
    esac
done

if [[ -z "$genome" ]]; then
    echo "Error: -g must be provided."
    exit 1
fi

mkdir -p ${output_dir}
cd ${output_dir} || exit 1

if [ ! -f "$lib_file" ]; then
    echo "Error: Library file not found: $lib_file"
    echo "Ensure you have run the full MCHelper automatic pipeline, or provide -l <library.fa>"
    exit 1
fi

# Run MCHelper
conda activate MCHelper
python3 $(which MCHelper.py) \
    -r T \
    --input_type fasta \
    -l "$lib_file" \
    -g "$genome" \
    -o . \
    -t "$THREADS" \
    -v Y > ${log_dir}/mchelper_manual.log
