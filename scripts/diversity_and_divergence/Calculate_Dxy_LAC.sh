#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --array=0-7
#SBATCH --output=Calculate_Dxy_LAC_%A_%a.oe
#SBATCH --job-name="CtDxyLAC"

source activate pixy

cd ../../processed_data/Dxy_LAC

VCF_DIR=$PWD/vcf
POP_DIR=$PWD
OUT_DIR=$POP_DIR/results
mkdir -p $OUT_DIR

sample_lists=($(ls -1 LAC_*.txt | sort))
POP_FILE=${sample_lists[$SLURM_ARRAY_TASK_ID]}

prefix=${POP_FILE%.txt}

VCFI=${VCF_DIR}/WI.20250627.hard-filter.isotype_${prefix}.vcf.gz

echo "Processing ${POP_FILE}"
echo "VCF: ${VCFI}"

pixy --stats dxy \
  --vcf $VCFI \
  --populations $POP_FILE \
  --output_folder $OUT_DIR \
  --window_size 10000 \
  --bypass_invariant_check yes \
  --output_prefix ${prefix}_ \
  --n_cores 15
