#!/bin/bash
# -----------------------------Name of the job-------------------------
#SBATCH --job-name=batch_mchelper
#-----------------------------Output files-----------------------------
#SBATCH --output=batch_mchelper.%j.out
#SBATCH --error=batch_mchelper.%j.err
#-----------------------------Required resources-----------------------
#SBATCH --time=36:00:00
#SBATCH -c 8
#SBATCH --mem 48GB
#-----------------------------Modules----------------------------------
module load miniconda3/4.9.2
eval "$(conda shell.bash hook)"

#which python

# Parse command-line options
while getopts ":s:n:l:g:b:h" option; do
    case "${option}" in
        s) species=${OPTARG} ;;
        n) strain=${OPTARG} ;;
        l) library=${OPTARG} ;;
        g) genome=${OPTARG} ;;
        b) BUSCO=${OPTARG} ;;
        h) echo "Usage: sbatch script.sh -s <species> -n <strain> -l <library.fa> -g <genome.fa> -b <BUSCO.hmm>"; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG"; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument."; exit 1 ;;
    esac
done

if [[ -z "$species" || -z "$strain" || -z "$library" || -z "$genome" || -z "$BUSCO" ]]; then
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
# By default MCHelper uses all available cores

conda activate MCHelper

python3 /mnt/netapp2/Store_csgcyjgp/dean/mctrials/TEammo/mchelper-ats/MCHelper.py \
-r A  \
-a F \
--input_type fasta \
-l $clean_lib \
-g $genome \
-b $BUSCO \
-o . \
-c 1 > mchelper_automatic.log

