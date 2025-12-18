#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Ct_run_admixture_post_processing.oe
#SBATCH --job-name="CtAdmPostPro"





cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data
mkdir -p Ct_admixture
cd Ct_admixture

out_folder="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Ct_admixture/"



echo "Starting post-processing..."
# merge log files
grep -h CV log_*.out | \
cut -f3- -d" " | \
sed 's/[(|):]//g' | \
sort -k1n | \
awk 'BEGIN{OFS="\t"; print "K", "CV"}; {print $0}' > admix_replicates_CV.tsv

echo "finish post-processing"



########################################################################
######### We no longer use this file to determine the best K ####

# # calculate avrage CV value
# source activate datamash
# datamash -g 1 mean 2 -H < "${out_folder}/admix_replicates_CV.tsv" | \
# sed 's/GroupBy(|)//g' > "${out_folder}/admix_summarized_CV.txt"

# # find best K value
# awk 'NR>1{print $0, sprintf("%3.2f", $2-p)} {p = $2}' "${out_folder}/admix_summarized_CV.txt" | \
# sed 's/-//g' | \
# awk '$3 == 0.00 {print $1-2, $1-1, $1, $1+1, $1+2}' > "${out_folder}/bestK.txt"

# echo "Post-processing completed. Best K candidates saved to bestK.txt"

######### We no longer use this file to determine the best K ####
########################################################################

