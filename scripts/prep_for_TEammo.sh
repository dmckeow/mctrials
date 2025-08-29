#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <Name for the data folder> <Sample csv file - with columns: species strain genome_path rm2_library_path curation2_final_lib busco_lib curation1 curation2 - see species_strains.csv>"
    exit 1
fi

data=$1
csv_file=$2


mkdir -p ${data} ${data}/RM2_output ${data}/MCH_output ${data}/0_raw ${data}/BUSCO_libs

# Prepare sample folders
# Loop through all subdirectories in data folder, excluding BUSCO_libs and db
for dir in ${data}/*/; do
    [[ "$(basename "$dir")" == "BUSCO_libs" || "$(basename "$dir")" == "db" ]] && continue

    # mkdir for species, strain
    while IFS=, read -r species strain genome_path rm2_library_path curation2_final_lib busco_lib curation1 curation2; do
        # Skip empty lines or comments
        [[ -z "$species" || "$species" == \#* ]] && continue
        mkdir -p "${dir}${species}/${strain}"
    done < "$csv_file"
done

# copy the genomes and libraries over
while IFS=, read -r species strain genome_path rm2_library_path curation2_final_lib busco_lib curation1 curation2; do
    # Skip empty lines or comments
    [[ -z "$species" || "$species" == \#* ]] && continue
    cp -n $genome_path ${data}/0_raw/${species}/${strain}
    cp -n $rm2_library_path ${data}/RM2_output/${species}/${strain}
done < "$csv_file"

# Get the BUSCO DB(s)
while IFS=, read -r species strain genome_path rm2_library_path curation2_final_lib busco_lib curation1 curation2; do
    # Skip empty lines or comments
    [[ -z "$species" || "$species" == \#* ]] && continue

    merged_hmm=${data}/BUSCO_libs/${busco_lib}_ALL.hmm
    if [ ! -f "$merged_hmm" ]; then
        echo "Getting busco lib ${busco_lib}"
        busco --download ${busco_lib}
        mv busco_downloads/lineages/${busco_lib}/ ${data}/BUSCO_libs/
        rm -fr busco_downloads/
        # The HMM profiles must be merged for MCHelper
        cat ${data}/BUSCO_libs/${busco_lib}/hmms/*.hmm > ${merged_hmm}
        rm -fr ${data}/BUSCO_libs/${busco_lib}
    else
        echo "BUSCO libs already exist, skipping"
    fi
done < "$csv_file"

# Prep commands to run MCHelper on them
# Make a plain file of the commands to submit for every genome and library in species_strains (for mctrials data):
while IFS=, read -r species strain genome_path rm2_library_path curation2_final_lib busco_lib curation1 curation2; do
    [[ -z "$species" || "$species" == \#* ]] && continue
    echo "cd ${data}/MCH_output/${species}/${strain}; \
        $(realpath scripts/mchelper.sh) \
        -s $species \
        -n $strain \
        -g $(realpath ${data}/0_raw/${species}/${strain}/$(basename $genome_path)) \
        -l $(realpath ${data}/RM2_output/${species}/${strain}/$(basename $rm2_library_path)) \
        -b $(realpath ${data}/BUSCO_libs/${busco_lib}_ALL.hmm); \
        cd -"
done < "$csv_file" | sed -E 's/ +/ /g' > batch_run_mchelper_${data}.sh


# make batch running file for this:
while IFS=, read -r species strain genome_path rm2_library_path curation2_final_lib busco_lib curation1 curation2; do
    [[ -z "$species" || "$species" == \#* ]] && continue
    echo "cd ${data}/MCH_output/${species}/${strain}; \
        $(realpath scripts/mchelper_teaid.sh) \
        -g $(realpath ${data}/0_raw/${species}/${strain}/$(basename $genome_path)); \
        cd -"
done < "$csv_file" | sed -E 's/ +/ /g' > batch_run_mchelper_teaid_${data}.sh