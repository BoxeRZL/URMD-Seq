#!/bin/bash

# === Auto-merge function ===

sample_file=$(ls -1 *.worker_*.csv | head -n 1)

if [[ -z "$sample_file" ]]; then
    echo "No worker output files found to merge."
    return 1
fi

outPathBase=${sample_file%.worker_*}

if [[ -f "$outPathBase" ]]; then
    mv "$outPathBase" "${outPathBase}.bak_$(date +%s)"
    echo "Existing $outPathBase backed up."
fi

header="groupSize,mutationFreq,lane,position,mid,reference,coverage,a,t,c,g,n,-,nonWTbases"
echo "$header" > "$outPathBase"

cat *.worker_*.csv >> "$outPathBase"

#rm -f *.worker_*.csv

echo "Merged all worker outputs into: $outPathBase"


