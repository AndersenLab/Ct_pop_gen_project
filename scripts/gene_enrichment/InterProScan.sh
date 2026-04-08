#!/bin/bash

#SBATCH -J IntProScan                 # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                      # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 48                            # Number of cores

wkdir="../../processed_data/gene_enrichment"

./interproscan.sh \
	--formats TSV \
	--input $wkdir/processed_data/gene_enrichment/c_tropicalis.NIC58_20251002.csq.longest.PC.noMtDNA.fa \
	--goterms \
	--cpu 48 \
	--iprlookup \
	--disable-precalc \
	--output-file-base $wkdir/NIC58_IPR_allApps_20251202 \
	--tempdir $TMP
