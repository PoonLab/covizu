#bash script to backup database, run ChunkyBot & update local flat file

#Copy database
cp data/gsaid.db "data/backup/$(date +"%Y_%m_%d").gsaid.db"

#Run updater script
python3 scripts/ChunkyBot.py -db "data/backup/$(date +"%Y_%m_%d").gsaid.db" >> debug/Chunkybot.log

#Migrate missing entries
python3 scripts/db_utils.py --targetdb "data/backup/$(date +"%Y_%m_%d").gsaid.db" >> debug/Migrate.log

#Delete old database files
cd data/backup
rm -v !("$(date +"%Y_%m_%d").gsaid.db")

# screen for non-human and low-coverage samples -> gisaid-filtered.fa
python3 scripts/filtering.py

# Run pangolin lineage assingment & upload to database
pangolin data/gisaid-filtered.fa -o data --outfile data/gisaid-filtered.pango.csv
python3 scripts/db_utils --lineagecsv data/gisaid-filtered.pango.csv

# calculate TN93 distances
tn93 -t 0.00005 -o data/gisaid.tn93.csv data/gisaid-filtered.fa

# cluster genomes into variants -> variants.csv, variants.fa
python3 scripts/variants.py

# calculate TN93 distances for clusters and output as HyPhy matrix
tn93 -o data/variants.tn93.txt -f hyphy data/variants.fa

# convert HyPhy matrix format to CSV
sed -i 's/[{}]//g' data/variants.tn93.txt

# hierarchical clustering -> data/clusters.json
Rscript scripts/hclust.R

# run FastTree and TreeTime
python3 scripts/treetime.py
