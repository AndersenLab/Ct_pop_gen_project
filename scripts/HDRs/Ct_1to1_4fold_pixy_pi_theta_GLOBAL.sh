#!/bin/bash
#SBATCH -A mschatz1_bigmem
#SBATCH -p bigmem
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=Ct_1to1_4fold_pixy_pi_theta_GLOBAL.oe
#SBATCH --job-name="Ct4foldPixyGLOBAL"

source activate pixy2026

cd ../../processed_data/HDR_stats

mkdir -p pi_theta_d_1to1_4fold_pixy

VCF="../../data/VCF/all.vcf.gz"
POP_file="../../processed_data/pi_theta_d/pop_file_GLOBAL.txt"

POP_DIR=$PWD
SITES_FILE="$POP_DIR/Ct_1to1_ortholog_4fold.1to1_4fold_sites.pixy_sites.txt"
OUT_DIR="$POP_DIR/pi_theta_d_1to1_4fold_pixy"

PREFIX="Ct_1to1_4fold_GLOBAL"

pixy \
  --stats pi watterson_theta \
  --vcf "$VCF" \
  --populations "$POP_file" \
  --sites_file "$SITES_FILE" \
  --window_size 10000 \
  --n_cores 24 \
  --output_folder "$OUT_DIR" \
  --output_prefix "$PREFIX"

echo "Done."
ls -lh pi_theta_d_1to1_4fold_pixy
