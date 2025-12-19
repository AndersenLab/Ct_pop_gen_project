#filex == file extension
#sp == species

filex="tsv"
sp="CT"
pattern='\_transformed\.tsv'

awk -v OFS='\t' '{print $0, FILENAME}' *.$filex | grep -v "^\[" | sed "s/${pattern}//" > ${sp}_all_nucmer_aln.tsv
