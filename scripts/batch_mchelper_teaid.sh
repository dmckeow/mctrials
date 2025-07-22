#!/bin/bash
# -----------------------------Name of the job-------------------------
#SBATCH --job-name=batch_mchelper_teaid
#-----------------------------Output files-----------------------------
#SBATCH --output=batch_mchelper_teaid.%j.out
#SBATCH --error=batch_mchelper_teaid.%j.err
#-----------------------------Required resources-----------------------
#SBATCH --time=36:00:00
#SBATCH -c 8
#SBATCH --mem 48GB
#-----------------------------Modules----------------------------------
module load miniconda3/4.9.2
eval "$(conda shell.bash hook)"

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

####
# Ok so we just have to TEAid independently through MCHelper - but only the version forked by Adrian
# Input must be the autocurated lib from MCHelper!

conda activate MCHelper
python3 /mnt/netapp2/Store_csgcyjgp/dean/mctrials/TEammo/mchelper-ats/MCHelper.py \
-r T \
--input_type fasta \
-l ../curated_sequences_NR.fa \
-g $genome \
-o . > ../mchelper_manual.log
