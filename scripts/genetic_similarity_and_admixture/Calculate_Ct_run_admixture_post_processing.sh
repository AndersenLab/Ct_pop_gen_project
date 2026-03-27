#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Ct_run_admixture_post_processing.oe
#SBATCH --job-name="CtAdmPostPro"


cd ../../processed_data
mkdir -p Ct_admixture
cd Ct_admixture

out_folder="./"

echo "Starting post-processing..."
# merge log files
grep -h CV log_*.out | \
cut -f3- -d" " | \
sed 's/[(|):]//g' | \
sort -k1n | \
awk 'BEGIN{OFS="\t"; print "K", "CV"}; {print $0}' > admix_replicates_CV.tsv
echo "finish post-processing"
