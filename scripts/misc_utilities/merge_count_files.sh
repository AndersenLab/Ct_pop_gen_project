awk -v OFS='\t' '{print $0, FILENAME}' *_counts.tsv | grep -v "^#" | sed 's/\.variant_counts.tsv//' > CT_all_strain_vc.tsv
