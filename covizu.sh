#!/bin/bash

# TODO: add download and updater scripts -> data/gisaid-aligned.fa

# screen for non-human and low-coverage samples -> gisaid-filtered.fa
python3 scripts/filtering.py

# calculate TN93 distances
tn93 -t 0.0001 -o data/gisaid.tn93.csv data/gisaid-filtered.fa

# cluster genomes into variants -> variants.csv, variants.fa
python3 scripts/variants.py

# calculate TN93 distances for clusters and output as HyPhy matrix
tn93 -o data/variants.tn93.txt -f hyphy data/variants.fa

# convert HyPhy matrix format to CSV
sed -i 's/[{}]//g' data/variants.tn93.txt

# hierarchical clustering -> data/clusters.json
Rscript scripts/hclust.R

# run FastTree and TreeTime
if [ ! -d "treetime" ]; then
  mkdir treetime
fi
python3 scripts/treetime.py
