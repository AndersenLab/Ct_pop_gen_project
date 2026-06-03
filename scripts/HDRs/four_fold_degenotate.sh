#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --output=FoldDegenotate.oe
#SBATCH --job-name="Ct4fold"

source activate degenotate

cd ./processed_data/HDR_stats/

GFF="../gene_enrichment/c_tropicalis.NIC58_20251002.csq.longest.gff3"
GENOME="../../data/c_tropicalis.NIC58_nanopore.June2021.genome.fa"

OUTDIR="degenotate_NIC58"
mkdir -p ${OUTDIR}

echo "GFF: ${GFF}"
echo "GENOME: ${GENOME}"
echo "OUTDIR: ${OUTDIR}"

ls -lh ${GFF}
ls -lh ${GENOME}

degenotate.py \
  -a ${GFF} \
  -g ${GENOME} \
  -o ${OUTDIR} \
  --overwrite

echo "Degenotate finished."
echo "Output files:"
find ${OUTDIR} -maxdepth 2 -type f -print
