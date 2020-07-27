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
