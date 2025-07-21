#!/bin/bash
# -----------------------------Name of the job-------------------------
#SBATCH --job-name=prep_ext_MCH
#-----------------------------Output files-----------------------------
#SBATCH --output=prep_ext_MCH.%j.out
#SBATCH --error=prep_ext_MCH.%j.err
#-----------------------------Required resources-----------------------
#SBATCH --time=24:00:00
#SBATCH -c 8
#SBATCH --mem 48GB
#-----------------------------Modules----------------------------------
module load miniconda3/4.9.2
eval "$(conda shell.bash hook)"

# Needs seqkit - i put it in base env
# Made to just run locally on WS1 in my home folder?
# Run in directory with MCHelper outputs run outsdie of TEAmmo

# Parse command-line options
while getopts ":l:g:s:h" option; do
    case "${option}" in
        l) library=${OPTARG} ;;
        g) genome=${OPTARG} ;;
        s) strain=${OPTARG} ;;        # To name clean families file
        h) echo "Usage: bash script.sh -l <library.fa> -g <genome.fa> -s <strain_name>"; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG"; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument."; exit 1 ;;
    esac
done

if [[ -z "$library" || -z "$genome" || -z "$strain" ]]; then
    echo "Error: All options -l, -g must be provided."
    exit 1
fi

# Prepare the clean library
clean_lib=$(echo $PWD"/"$strain"-clean_families.fa")
conda activate bioinf # This is where i keep common tools
awk '/^>/' $library | \
awk '!/Satellite|Simple_repeat|tRNA|rRNA|Retroposon|snRNA|scRNA/' $library | \
sed 's/^>//g' | \
seqkit grep -n -f - $library > $clean_lib

#log=$(ls -t *.out | head -1) 
#cp $log mchelper_automatic.log # made by redirect in first script instead

# Close, can see the listed files, it just freezes
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
