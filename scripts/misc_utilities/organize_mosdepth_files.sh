mkdir regions
mkdir misc
mkdir thresh

mv *.csi misc/
mv *.txt misc/
mv *.regions.* regions/
mv *.thresholds.* thresh/
cd thresh
gunzip *
awk -v OFS='\t' '{print $0, FILENAME}' *.bed | grep -v "^#" | sed 's/\.thresholds.bed//' > CT_all_thresh_cov.tsv
