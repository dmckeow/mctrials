#!/bin/bash
# -----------------------------Name of the job-------------------------
#SBATCH --job-name=batch_mchelper
#-----------------------------Output files-----------------------------
#SBATCH --output=batch_mchelper.%j.out
#SBATCH --error=batch_mchelper.%j.err
#-----------------------------Required resources-----------------------
#SBATCH --time=48:00:00
#SBATCH -c 8
#SBATCH --mem 64GB
#-----------------------------Modules----------------------------------
module load miniconda3/4.9.2
eval "$(conda shell.bash hook)"
conda activate MCHelper

#which python

# Parse command-line options
while getopts ":s:n:l:g:b:h" option; do
    case "${option}" in
        s) sp_name=${OPTARG} ;;       # -s for species name
        n) strain_name=${OPTARG} ;;   # -n for strain name
        l) library=${OPTARG} ;;       # -l for library file
        g) genome=${OPTARG} ;;        # -g for genome file
        b) BUSCO=${OPTARG} ;;         # -b for BUSCO file
        h) echo "Usage: sbatch script.sh -s <species> -n <strain> -l <library.fa> -g <genome.fa> -b <BUSCO.hmm>"; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG"; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument."; exit 1 ;;
    esac
done

if [[ -z "$sp_name" || -z "$strain_name" || -z "$library" || -z "$genome" || -z "$BUSCO" ]]; then
    echo "Error: All options -s, -n, -l, -g, and -b must be provided."
    exit 1
fi

# Prepare the clean library
clean_lib=$(echo $PWD"/"$strain"-clean_families.fa")
conda activate bioinf # This is where i keep common tools
awk '/^>/' $library | \
awk '!/Satellite|Simple_repeat|tRNA|rRNA|Retroposon|snRNA|scRNA/' $library | \
sed 's/^>//g' | \
seqkit grep -n -f - $library > $clean_lib

# The outputs go to pwd
#source activate MCHelper
# By default MCHelper uses all available cores

#python /mnt/netapp2/Store_csgcyjgp/dean/mctrials/MCHelper/MCHelper.py \
python3 /mnt/netapp2/Store_csgcyjgp/dean/mctrials/TEammo/mchelper-ats/MCHelper.py \
-r A  \
-a F \
--input_type fasta  \
-l $clean_lib  \
#-l $library  \
-g $genome \
-b $BUSCO  \
-o . \
-c 1 