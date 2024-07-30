#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=download_and_unzip_tif_files.oe
#SBATCH --job-name=dl_uz_tif_files


cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/
mkdir -p tif
cd tif

wget https://figshare.com/ndownloader/articles/16571064/versions/6


for year in {2000..2020}; do
    unzip "hfp${year}.zip"
done
