#!/bin/bash
#SBATCH -A mschatz1_bigmem
#SBATCH -p bigmem
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=Calculate_Dxy_LAC.oe
#SBATCH --job-name="CtDxyLAC"

source activate pixy2026

cd ../../processed_data/
mkdir -p Dxy_LAC
cd Dxy_LAC

vcf_raw="../../data/VCF/all.vcf.gz"

POP_DIR=$PWD
OUT_DIR=$POP_DIR/results
mkdir -p $OUT_DIR
PREFIX="Ct_Dxy_"

pixy \
 --stats pi dxy fst watterson_theta tajima_d \
 --vcf "$vcf_raw" \
 --populations "$POP_DIR/pop_file_lineage.txt" \
 --window_size 10000 \
 --n_cores 48 \
 --output_folder "$OUT_DIR" \
 --output_prefix "$PREFIX"

