# mctrials
Practice for manual curation of TEs in Drosophila
## Workflow

### Setup
Working temporarily in /home/csic/gcy/jgp/extra_storage/dean/mctrials/mctrials

```{bash}
conda create -n bioinf -y
conda activate bioinf
conda install -c bioconda seqkit
```

```{bash}
cd /home/csic/gcy/jgp/extra_storage/dean/mctrials/mctrials
mkdir -p data data/RM2_output data/MCH_output data/0_raw data/BUSCO_libs

# Strain folders
# Path to your CSV file
csv_file="species_strains.csv"
sed -i -z 's/$/\n/1' $csv_file

# Loop through all subdirectories in data/, excluding BUSCO_libs
for dir in data/*/; do
    [[ "$(basename "$dir")" == "BUSCO_libs" ]] && continue

    # mkdir for species, strain
    while IFS=, read -r species strain genome_path rm2_library_path busco_lib; do
        # Skip empty lines or comments
        [[ -z "$species" || "$species" == \#* ]] && continue
        mkdir -p "${dir}${species}/${strain}"
    done < "$csv_file"
done

# copy the genomes and libraries over
while IFS=, read -r species strain genome_path rm2_library_path busco_lib; do
    # Skip empty lines or comments
    [[ -z "$species" || "$species" == \#* ]] && continue
    cp -n $genome_path data/0_raw/${species}/${strain}
    cp -n $rm2_library_path data/RM2_output/${species}/${strain}
done < "$csv_file"

# Get the BUSCO DB
busco --download diptera_odb12
mv busco_downloads/lineages/diptera_odb12/ data/BUSCO_libs/
rm -fr busco_downloads/
# The HMM profiles must be merged for MCHelper
cat data/BUSCO_libs/diptera_odb12/hmms/*.hmm > data/BUSCO_libs/diptera_odb12_ALL.hmm
rm -fr data/BUSCO_libs/diptera_odb12


```

Then put genome and RM2 library in corresponding folder in RM2

### Running MCHelper
As batch jobs on CESGA, clunky maybe I can make a Nextflow module for this later...
* D.merina
* D.miranda
* D.santomea
* D.tristis
* C.elegans (toy test data)

```{bash}

module load miniconda3/4.9.2 # CESGA requirement

# Manual version with toy data:
cd data/MCH_output/C.elegans/N2_sub3
sbatch ../../../../scripts/mchelper.sh \
-s C.elegans \
-n N2_sub3 \
-l ../../../RM2_output/C.elegans/N2_sub3/N2_sub3-families.fa \
-g ../../../0_raw/C.elegans/N2_sub3/N2_subset3.fna \
-b ../../../BUSCO_libs/nematoda_odb10_ALL.hmm
cd -

# Make a plain file of the commands to submit for every genome and library in species_strains (for mctrials data):
while IFS=, read -r species strain genome_path rm2_library_path busco_lib; do
    [[ -z "$species" || "$species" == \#* ]] && continue
    echo "cd data/MCH_output/${species}/${strain}; \
        sbatch $(realpath scripts/mchelper.sh) \
        -s $species \
        -n $strain \
        -g $(realpath data/0_raw/${species}/${strain}/$(basename $genome_path)) \
        -l $(realpath data/RM2_output/${species}/${strain}/$(basename $rm2_library_path)) \
        -b $(realpath data/BUSCO_libs/${busco_lib}_ALL.hmm); \
        cd -"
done < "$csv_file" | sed -E 's/ +/ /g' > batch_run_mchelper_sh

```



### Running TEAid through MCHelper
* Need to run TEAid independently through MCHelper - but only the version forked by Adrian
* The library input MUST be the NR curated library output from MCHelper

```{bash}

cd /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/data/MCH_output/C.elegans/N2_sub3
sbatch /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/scripts/mchelper_teaid.sh \
-g /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/data/0_raw/C.elegans/N2_sub3/N2_subset3.fna

cd /home/csic/gcy/jgp/extra_storage/dean/mctrials/mctrials/data/MCH_output/D.miranda/v2.1_MSH22_RefSeq
sbatch /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/scripts/mchelper_teaid.sh \
-g /home/csic/gcy/jgp/extra_storage/dean/mctrials/mctrials/data/0_raw/D.miranda/v2.1_MSH22_RefSeq/D.miranda_v2.1_MSH22_RefSeq.fasta

cd /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/data/MCH_output/D.santomea/STO_CAGO_1482_RefSeq
sbatch /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/scripts/mchelper_teaid.sh \
-g /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/data/0_raw/D.santomea/STO_CAGO_1482_RefSeq/D.santomea_STO_CAGO_1482_RefSeq.fasta

cd /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/data/MCH_output/D.tristis/nanopore_D2
sbatch /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/scripts/mchelper_teaid.sh \
-g /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/data/0_raw/D.tristis/nanopore_D2/D.tristis_nanopore_D2.fasta

cd /home/csic/gcy/jgp/extra_storage/dean/mctrials/mctrials/data/MCH_output/D.merina/NA
sbatch /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/scripts/mchelper_teaid.sh \
-g /mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/data/0_raw/D.merina/NA/D.merina.rm.fasta

```

### Getting data into TEammo
I ran MCHelper on pregenerated raw RM2 libraries, using the CESGA cluster - now I will transfer the results to a local workstation

```{bash}
rsync -avz csgcyjgp@ft3.cesga.es:/mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/data/ mctrials_data/
```

Then we linked the data to the data folder where TEammo looks and stores outputs
* The toy data is linked with `ln -s ../teammo_data/ data`
```{bash}
cd ~/TEammo
rm -f data
ln -s ../mctrials_data/ data
```


### Running TEammo
* On a local workstation, ssh through own computer
* Login and installation setup by Adrian, but likely it was the conda setup for TEammo he did
  
```{bash}
ssh -L 3001:localhost:3001 dmckeow@161.111.135.77

cd TEammo
# If you just run the Shiny app in terminal without redirect it is unstable
Rscript TEammo_app.R > teammo.log 2>&1 &

# In case a process was not cancelled properly and it is stilling running
lsof -i :3001
kill -9 <pid>
```