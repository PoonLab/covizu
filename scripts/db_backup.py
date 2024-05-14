import subprocess
from datetime import datetime
import os
import glob
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="CoVizu backend database backup script")

    parser.add_argument('--outdir', type=str, default="/data/covizu_db_backups",
                        help='Directory to store the database backups')
    parser.add_argument('--dbname', type=str, default=os.environ.get("POSTGRES_DB", "gisaid_db"),
                        help="Postgresql database name")
    parser.add_argument('--dbuser', type=str, default=os.environ.get("POSTGRES_USER", None),
                        help="Postgresl user")
    parser.add_argument('--docker', type=str, default=os.environ.get("DOCKER_CONTAINER", "postgres"),
                        help="Docker container name")
    parser.add_argument('--data-dump', dest='datadump', action="store_true",
                        help="Perform a backup of the database")
    parser.add_argument('--remove-backup', dest='rmbkup', action="store_true",
                        help="Remove old database backups")

    return parser.parse_args()


def get_week_number(year, month, day):
    dt = datetime(year, month, day)
    first_day_of_month = datetime(dt.year, dt.month, 1)
    offset = (first_day_of_month.weekday() + 1) % 7  
    return (dt.day + offset - 1) // 7 + 1


def dump_psql_data(db_name, username, file_path, docker_container_name):
    current_date = datetime.now()
    year = current_date.year
    month = current_date.month
    day = current_date.day
    
    week_number = get_week_number(year, month, day)
    
    backup_file_name = f"db_backup.{year}-{month}-{day}-{week_number}.db"
    backup_file_path = os.path.join(file_path, backup_file_name)

    try:
        dump_command = f"docker exec -it {docker_container_name} pg_dump -U {username} -d {db_name} > {backup_file_path}"
        process = subprocess.Popen(dump_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        _, error = process.communicate()

        if process.returncode == 0:
            print(f"Data dumped successfully to {backup_file_path}")
        else:
            print(f"Error dumping data: {error.decode('utf-8')}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")


def remove_backups(directory):
    backups = glob.glob(os.path.join(directory, '*.db*'))
    if len(backups) == 0:
        print("No backups found")
        return

    for file in backups:
        file_name = os.path.basename(file)
        year, month, day, week_number = file_name.split('.')[1].split('-')

        if (datetime.now() - datetime(int(year), int(month), int(day))).days < 60 or int(week_number) == 1:
            continue

        try:
            os.remove(file)
            print(f"Removing {file}")
        except FileNotFoundError:
            print(f"Database backup already removed: {file}")


if __name__ == "__main__":
    args = parse_args()
    args.outdir = os.path.normpath(args.outdir)

    if args.datadump:
        dump_psql_data(args.dbname, args.dbuser, args.outdir, args.docker)

    if args.rmbkup:
        remove_backups(args.outdir)