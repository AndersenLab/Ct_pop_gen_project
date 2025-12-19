faidx="c_tropicalis.NIC58_nanopore.June2021.genome.fa.fai"
outfile="ONT_NIC58_1kb_bins.bed"

bedtools makewindows -g $faidx -w 1000 > $outfile 
