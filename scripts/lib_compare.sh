#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 8 ]; then
    echo "Usage: $0 lib1 lib1_name lib2 lib2_name species strain RM2_lib <MCHelper auto curated lib>"
    exit 1
fi

lib1=$1       # e.g. dean.fa
lib1_name=$2
lib2=$3       # e.g. marta.fa
lib2_name=$4
species=$5    # e.g. D.tristis
strain=$6     # e.g. nanopore_D2
RM2_lib=$7    # user-provided label
MCH_lib=$8    # user-provided label


outdir="results/lib_compare"
mkdir -p "${outdir}"


# Make BLAST databases
makeblastdb -in "${lib1}" -dbtype "nucl"
makeblastdb -in "${lib2}" -dbtype "nucl"
makeblastdb -in "${RM2_lib}" -dbtype "nucl"
makeblastdb -in "${MCH_lib}" -dbtype "nucl"

# Reciprocal BLAST
blastn -query "${lib1}" -db "${lib2}" -outfmt 6 -max_hsps 1 \
    -out "${outdir}/${species}_${strain}_${lib1_name}_vs_${lib2_name}.blast.out"

blastn -query "${lib2}" -db "${lib1}" -outfmt 6 -max_hsps 1 \
    -out "${outdir}/${species}_${strain}_${lib2_name}_vs_${lib1_name}.blast.out"

# BLAST lib2 against RM2 and MCH to recover original names
blastn -query "${lib2}" -db "${RM2_lib}" -outfmt 6 -max_hsps 1 \
    -out "${outdir}/${species}_${strain}_${lib2_name}_vs_RM2.blast.out"

blastn -query "${lib2}" -db "${MCH_lib}" -outfmt 6 -max_hsps 1 \
    -out "${outdir}/${species}_${strain}_${lib2_name}_vs_MCHelperAuto.blast.out"
