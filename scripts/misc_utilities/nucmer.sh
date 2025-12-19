#!/bin/bash

#activate environment
source activate nucmer

# align with nucmer (will spit out a .delta file)
nucmer --maxgap=500 --mincluster=100 --prefix=$prefix --coords $ref $query

# get coordinate file
show-coords -r -l -T $prefix.delta | awk '$5 > 1000' > ${prefix}_transformed.tsv 
