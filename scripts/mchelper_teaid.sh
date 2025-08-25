#!/bin/bash
eval "$(conda shell.bash hook)"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
echo "Using ${THREADS} cores"


# Needs seqkit - i put it in base env
# Made to just run locally on WS1 in my home folder?
# Run in directory with MCHelper outputs run outsdie of TEAmmo

# Parse command-line options
while getopts ":g:h" option; do
    case "${option}" in
        g) genome=${OPTARG} ;;
        h) echo "Usage: bash script.sh -g <genome.fa>"; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG"; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument."; exit 1 ;;
    esac
done

if [[ -z "$genome" ]]; then
    echo "Error: -g must be provided."
    exit 1
fi

mkdir -p 1_MI_MCH
cd 1_MI_MCH

if [ ! -f ../curated_sequences_NR.fa ]; then
    echo "No auto curated library from MCHelper found in parent directory:"
    echo "Ensure you have run the full MCHelper automatic pipeline, and are running this from its output folder"
fi

####
# Ok so we just have to TEAid independently through MCHelper - but only the version forked by Adrian
# Input must be the autocurated lib from MCHelper!

conda activate MCHelper
python3 ~/tools/TEammo/mchelper-ats/MCHelper.py \
-r T \
--input_type fasta \
-l ../curated_sequences_NR.fa \
-g $genome \
-o . \
-t $THREADS \
-v Y > ../mchelper_manual.log

# -v Y added after the library size discrepancy noticed