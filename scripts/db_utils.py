import sqlite3
import argparse
import csv
import sys
from tempfile import NamedTemporaryFile


def export_fasta(cursor, outfile='data/gisaid-aligned.fa'):
    """ Function that writes fasta to outfile
    :params:
        :cursor: sqlite3 object
        :outfile: (str) output file
    """
    data = cursor.execute('SELECT `header`, `aligned` FROM SEQUENCES;').fetchall()
    out = open(outfile, 'w+')
    for h, s in data:
        out.write('>{}\n{}\n'.format(h, s))
    out.close()


def open_connection(database):
    """ open connection to database, initialize tables if they don't exist
        :params:
            :database: str, name of database
        :out:
            :cursor: interactive sql object containing tables
    """
    conn = sqlite3.connect(database, check_same_thread=False)
    cur = conn.cursor()

    # create tables if they don't exist
    seqs_table = 'CREATE TABLE IF NOT EXISTS SEQUENCES (accession VARCHAR(255) ' \
                 'PRIMARY KEY, header VARCHAR(255), aligned BLOB);'
    cur.execute(seqs_table)

    sample_table = "CREATE TABLE IF NOT EXISTS SAMPLE (accession VARCHAR(255) " \
                   "PRIMARY KEY, virus_name VARCHAR(255), collection_date VARCHAR(10), " \
                   "location VARCHAR(30), lineage VARCHAR(20), clade VARCHAR(30) );"
    cur.execute(sample_table)

    sequencing_table = "CREATE TABLE IF NOT EXISTS SEQUENCING (accession VARCHAR(255) " \
                       "PRIMARY KEY, platform VARCHAR(100), assembly_method VARCHAR(100), " \
                       "sample_type VARCHAR(100))"

    cur.execute(sequencing_table)

    acknowledgement_table = "CREATE TABLE IF NOT EXISTS ACKNOWLEDGEMENTS (accession " \
                            "VARCHAR(255) PRIMARY KEY, originating_lab VARCHAR(100), " \
                            "sub_lab VARCHAR(100), authors BLOB)"
    cur.execute(acknowledgement_table)

    ptinfo_table = "CREATE TABLE IF NOT EXISTS PTINFO (accession VARCHAR(255) PRIMARY KEY, " \
                   "gender VARCHAR(30), age INT, ptstatus VARCHAR(20))"
    cur.execute(ptinfo_table)

    raw_seq_table = "CREATE TABLE IF NOT EXISTS RAWSEQ (accession VARCHAR(255) PRIMARY KEY, " \
                    "header VARCHAR(255), unaligned BLOB);"
    cur.execute(raw_seq_table)

    lineage_table = "CREATE TABLE IF NOT EXISTS LINEAGE (entry_id INTEGER PRIMARY KEY, " \
                    "accession VARCHAR(255), lineage varchar(15), probability REAL, pangoLEARN_version DATE," \
                    "status VARCHAR(255), note BLOB);"

    cur.execute(lineage_table)

    conn.commit()
    return cur, conn


def insert_seq(cursor, header, aligned):
    """ Wrapper function for aligning & inserting into SEQUENCES table
        :params:
            :cursor: sqlite database handler
            :seq: raw sequence
            :header: fasta header
            :aligned: aligned sequence
    """
    vars= [header.split('|')[1], header, aligned]
    result = cursor.execute(
        "REPLACE INTO SEQUENCES('accession', 'header', 'aligned') "
        "VALUES(?, ?, ?)", vars)
    return 0


def insert_into_rawseqs(database, fasta):
    """ Wrapper function for inserting into RAWSEQ table
        :params:
            :database: sqlite database 
            :fasta: file containing raw seqs
    """
    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    handle = convert_fasta(open(fasta))

    inserted_acessions = []
    for header,seq in handle:
        vars= [header.split('|')[1], header, seq]
        result = cursor.execute(
            "REPLACE INTO RAWSEQ ('accession', 'header', 'unaligned') VALUES "
            "(?,?,?)", vars)
        inserted_acessions.append(vars[0])
    conn.commit()
    conn.close()
    return inserted_acessions


def iterate_lineage_csv(cursor, csvFile):
    """
    Wrapper function for inserting csv's into LINEAGES table
    :params:
        :cursor: sqlite database handler
        :csvFile: output of pangolin containing run data
    :out:
        :None:

    """
    with open(csvFile, "r") as infile:
        reader = csv.reader(infile)
        next(reader, None) # skip header
        for row in reader:
            vars = [row[0].split('|')[1], row[1], row[2], row[3], row[4], row[5]]
            cursor.execute(
                "INSERT INTO LINEAGE(`accession`, `lineage`, `probability`, "
                "`pangoLEARN_version`, `status`, `note`) VALUES(?,?, ?, DATE(?),?,?)",
                vars
            )


def insert_lineage(cursor, payload):
    """
    Wrapper function for streaming Lineage typing into sqlite database -Lineage table-
        :params:
            :cursor: sqlite database handler
            :payload: data streamed from pangorider [accession, lineage, probability, pangoLEARN_version ]
        :output:
            None
    """
    cursor.execute(
        "INSERT INTO LINEAGE(`accession`, `lineage`, `probability`, "
        "`pangoLEARN_version`, `status`, `note`) VALUES(?,?, ?, DATE(?),?,?)",
        payload
    )


def process_tech_meta(cursor, tsvFile):
    """
    Wrapper function for inserting into SAMPLE & SEQUENCING table
    #FROM SEQUENCING TECH METADATA
    :params:
        :cursor: sqlite database handler
        :tsvFile: tabs separated values file containing data, from the
        Sequencing technology metadata download option
    :out:
        :debug: report
    """

    with open(tsvFile) as handle:
        reader= csv.reader(handle, delimiter="\t")
        header = next(reader)

        for row in reader:
            sample_vars = [row[1], row[0], row[2], row[3], row[12], row[13]]
            result1 = cursor.execute(
                "REPLACE INTO SAMPLE('accession', 'virus_name', "
                "'collection_date', 'location', 'lineage', 'clade') VALUES(?, ?, ?, ?, ?, ?)",
                sample_vars
            )
            sequencing_vars = [row[1], row[8], row[9], row[6]]
            result2 = cursor.execute(
                "REPLACE INTO SEQUENCING('accession', 'platform', 'assembly_method', "
                "'sample_type') VALUES(?, ?, ?, ?)",
                sequencing_vars
            )


def insert_acknowledgement(cursor, tsvFile):
    """
    Wrapper function for processing tsvFile and inserting into PTINFO table
    #FROM Patient status metadata
    :params:
        :cursor: sqlite database handler
        :tsvFile: tabs separated values file containing data, from the
        Sequencing technology metadata download option
    :out:
        :debug: report
    """

    with open(tsvFile) as handle:
        reader= csv.reader(handle, delimiter="\t")
        header = next(reader)

        for row in reader:
            vars = [row[0], row[4], row[5], row[6]]
            result = cursor.execute(
                "REPLACE INTO ACKNOWLEDGEMENTS('accession', 'originating_lab', 'sub_lab', "
                "'authors') VALUES(?, ?, ?, ?)",
                vars
            )
    handle.close()


def process_ptinfo(cursor, tsvFile):
    """
    Wrapper function for processing tsvFile and inserting into PTINFO table
    #FROM Patient status metadata
    :params:
        :cursor: sqlite database handler
        :tsvFile: tabs separated values file containing data, from the patient
        status metadata download option
    :out:
        :debug: report
    """

    with open(tsvFile) as handle:
        reader= csv.reader(handle, delimiter="\t")
        header = next(reader)
        for row in reader:
            #PTINFO table
            vars = [row[1], row[5], row[6], row[7]]
            result = cursor.execute(
                "REPLACE INTO PTINFO('accession', 'gender', 'age', 'ptstatus') "
                "VALUES(?, ?, ?, ?)",
                vars
            )
    handle.close()


def pull_field(cursor, field):
    """
    Wrapper function to return dictionary key-val pair being accession-field
    """
    field = ["`" + field + "`"]
    result = cursor.execute("SELECT `accession`, ? FROM SEQUENCES;", field)
    field_dict= dict(result.fetchall())
    return field_dict


def report_changes(database, header, seq):
    """
    Function to check for changes in uploaded seqs, replaces different headers
    :params:
        :cursor:
        :header: seq name
        :seq: unaligned covid seq
    """
    cursor, conn = open_connection(database) 

    accession = header.split('|')[1]
    record = cursor.execute('SELECT header, unaligned from rawseq WHERE accession = ?',
                            [accession]).fetchone()
    if type(record) == None:
        conn.commit()
        conn.close()
        return [header, seq]

    old_header, old_seq = record #otherwise unpack values

    if header != old_header:
        cursor.execute('REPLACE INTO sequences SET header = ? WHERE accession = ?',
                       [header, accession])
        conn.commit()
        conn.close()
        return 1
    if hash(seq) != hash(old_seq): #hash comparison is quicker?
        conn.commit()
        conn.close()
        return [header, seq]
    conn.commit()
    conn.close()

    return 0


def convert_fasta(handle):
    """
    Taken from gotoh2 by artpoon
    Parse FASTA file as a list of header, sequence list objects
    :param handle:  open file stream
    :return:  List of [header, sequence] records
    """
    result = []
    h, sequence = None, ''

    for line in handle:
        if line.startswith('$'):  # skip comment
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h, sequence])
                sequence = ''
            h = line.lstrip('>').rstrip()
        else:
            sequence += line.strip().upper()

    result.append([h, sequence])  # handle last entry
    return result


def iterate_fasta(fasta, database = 'data/gsaid.db'):
    """
    Function to iterate through aligned fasta file *non threaded*
    :param cursor: sqlite database handler
    :param fasta: file containing sequences
    :param ref: file containing reference sequence, default-> NC_04552.fa
    """
    handle = convert_fasta(open(fasta))
    cursor, conn = open_connection(database)

    for header, seq in handle:
        insert_seq(cursor, header, seq)
        conn.commit()
    conn.close()


def iterate_handle(handle, database='data/gsaid.db'):
    """
    Function to iterate through list [(header,seq)....]
    :param cursor: sqlite database handler
    :param handle: list containing tuples (header, raw seq)
    :param ref: file containing reference sequence, default-> NC_04552.fa
    """

    cursor, conn = open_connection(database)

    # align and insert into database
    for header, alignedseq in handle:
        insert_seq(cursor, header, alignedseq)
    conn.commit()


def migrate_entries(old_db, new_db):
    """
    Function to migrate missing seqs from old_db to new_db
    :params:
        :old_db: sqlite3 target db
        :new_db: sqlite3 db with missing seqs
    """
    oldconnect = sqlite3.connect(old_db)
    chunkyconnect = sqlite3.connect(new_db)
    oldcursor = oldconnect.cursor()
    chunkcursor = chunkyconnect.cursor()

    # get existing accessions in target db, create list to derive missing accesions
    existing_accessions = list(oldcursor.execute('SELECT accession from SEQUENCES;').fetchall())
    EA_list = []
    for accession in existing_accessions:
        EA_list.append(accession[0])

    # get all sequences from the source db and insert them into target
    all_chunk_seqs = list(chunkcursor.execute('SELECT * FROM SEQUENCES;').fetchall())
    count = 0
    for accession, header, unaligned_hash, aligned in all_chunk_seqs:
        if accession not in EA_list:
            count+=1
            oldcursor.execute('INSERT INTO SEQUENCES (accession, header, unaligned_hash, aligned) '
                              'VALUES (?,?,?,?);', [accession, header, unaligned_hash, aligned])
            oldconnect.commit()

    print('Updated {}, inserted {} seqs.'.format(old_db, str(count)))

    # get all existing hashes from target db, create list to derive missing accesions
    existing_hashes = list(oldcursor.execute('SELECT hash_key FROM HASHEDSEQS').fetchall())
    EH_list = []
    for hash_key in existing_hashes:
        EH_list.append(hash_key[0])

    hash_count = 0
    all_chunk_hashes = list(chunkcursor.execute('SELECT * FROM HASHEDSEQS;').fetchall())
    for hash_key, aligned_seq in all_chunk_hashes:
        if hash_key not in EH_list:
            hash_count+=1
            oldcursor.execute('INSERT INTO HASHEDSEQS (hash_key, aligned_seq) VALUES (?,?)',
                              [hash_key, aligned_seq])
            oldconnect.commit()

    print('Updated {}, inserted {} hashes.'.format(old_db, str(hash_count)))


def dump_raw(outfile, db='data/gsaid.db'):
    """
    Function to dump all raw seqs in RAQSEQ table
    :params:
        :outfile: destination (fasta) file to write seqs to
        :db: sqlite3 database
    """
    cur, conn = open_connection(db)
    seqs = cur.execute("SELECT `header`, `unaligned` FROM rawseq;").fetchall()
    with open(outfile, 'w') as fastafile:
        for h, s in seqs:
            fastafile.write('>{}\n{}\n'.format(h, s))
    fastafile.close()
    conn.close()


def retrieve_seqs(headers, db="data/gsaid.db"):
    """ Retrieve genomes with list of headers (labels) """
    cur, conn = open_connection(db)
    seqs = []
    for h in headers:
        result = cur.execute(
            "SELECT `aligned` from SEQUENCES where `header`='{}';".format(h)
        ).fetchall()
        # TODO: warn if multiple results?
        if len(result) == 0:
            print("ERROR: query {} failed to retrieve any genomes".format(headers))
            sys.exit()

        seqs.append(result[0][0])
    conn.close()
    return seqs


def dump_lineages(db='data/gsaid.db', by_lineage=False):
    """
    Connect to sqlite3 database and retrieve all records from the LINEAGE
    table.  Construct a dictionary keyed by accession number.
    :param db:  str, path to sqlite3 database
    :param by_lineage:  bool, if True then returned dict is keyed by lineage
    :return:  dict, lineage assignments keyed by accession (default), or lists of
              genome labels (headers) keyed by lineage if by_lineage is True
    """
    cursor, conn = open_connection(db)
    labels = cursor.execute("SELECT `accession`, `header` FROM SEQUENCES").fetchall()
    labels = dict(labels)

    data = cursor.execute("SELECT `accession`, `lineage`, `probability`, `pangoLEARN_version` "
                          "FROM LINEAGE;").fetchall()
    result = {}
    if by_lineage:
        # used by treetime.py
        for accn, lineage, prob, version in data:
            header = labels.get(accn, None)
            if header is None:
                continue
            if lineage not in result:
                result.update({lineage: []})
            result[lineage].append(labels[accn])
    else:
        for accn, lineage, prob, version in data:
            result.update({accn: {
                'lineage': lineage, 'prob': prob, 'version': version
            }})
    conn.close()
    return result


def dump_raw_by_lineage(db='data/gsaid.db'):
    """
    Creates a generator that yields handles to named temporary files
    containing the unaligned genomes for each Pangolin lineage
    :param db:  str, path to sqlite3 database
    :yield:  tuples, Pangolin lineage and path to temporary FASTA file
    """
    cursor, conn = open_connection(db)
    lineages = cursor.execute("select DISTINCT(lineage) from LINEAGE;").fetchall()
    for lineage in lineages:
        # retrieve unaligned genomes from db
        raw = cursor.execute(
            "select header, unaligned from RAWSEQ where accession in \
            (select accession from LINEAGE where lineage='{}');".format(lineage[0])
        ).fetchall()
        handle = NamedTemporaryFile('w', delete=False)
        for h, s in raw:
            handle.write('>{}\n{}\n'.format(h, s))
        handle.close()
        yield lineage[0], handle.name
    conn.close()


def process_meta(db, metafiles):
    """
    Function to process tsv files containing meta data
    :params:
        :db: sqlite3 database
        :metafiles: [meta_epi,meta_tech] - list containing absolute paths
    """
    cursor, conn = open_connection(db)
    # process pt_info file
    process_ptinfo(cursor, metafiles[0])

    # process tech file
    process_tech_meta(cursor, metafiles[1])
    conn.commit()
    conn.close()


def dump_concordance(db, concordanceOUT, discordanceOUT):
    """
    Function to dump lineage assignments with accession no's

    :params:
        :db: sqlite3 database
        :concordanceOUT: csv file for concordant records (same assignment gisaid/ local pangolin)
        :discordanceOUT: csv file for discordant records (diff assignment gisaid/ local pangolin)
    """
    cur, conn = open_connection(db)
    # Handle concordance first
    CONoutfile = open('data/{}.csv'.format(concordanceOUT), 'w')
    # take most recent lineage assigment, join lineage and sample tables on accession
    CONhandle = cursor.execute(
        'SELECT filtered_lineage.accession, filtered_lineage.lineage as local_lineage, probability, '
        'pangoLEARN_version, status, note, virus_name, collection_date, location, sample.lineage '
        'AS gisaid_lineage, clade '
        'FROM (SELECT accession, lineage, probability, max(pangoLEARN_version) '
        'AS pangoLEARN_version, status, note '
        'FROM LINEAGE GROUP BY accession) AS filtered_lineage, sample '
        'WHERE sample.accession = filtered_lineage.accession '
        'AND local_lineage == gisaid_lineage;'
    ).fetchall()

    # create start of csv file
    CONouthandle = ['accession,local_lineage,local_probability,local_pangoLEARN_version,local_status,seq_name,sample_col_date,sample_location,gisaid_lineage,gisaid_clade \n']
    # iterate handle, print out fields to csv
    for accession1, lineage1, prob, version, status, note, name, date, location, lineage2, clade in CONhandle:
        CONouthandle.append('{},{},{},{},{},{},{},{},{},{},{}\n'.format( accession1, lineage1, prob, version, status, note, name, date, location, lineage2, clade))
    CONoutfile.writelines(CONouthandle)  # process and write outfile
    CONoutfile.close()

    # repeat for dicordance
    DIShandle = cursor.execute(
        'SELECT filtered_lineage.accession, filtered_lineage.lineage '
        'as local_lineage, probability, pangoLEARN_version, status, note, virus_name, '
        'collection_date, location, sample.lineage '
        'AS gisaid_lineage, clade '
        'FROM (SELECT accession, lineage, probability, max(pangoLEARN_version) '
        'AS pangoLEARN_version, status, note '
        'FROM LINEAGE GROUP BY accession) AS filtered_lineage, sample '
        'WHERE sample.accession = filtered_lineage.accession AND local_lineage != gisaid_lineage;'
    ).fetchall()

    # open outfile, write header row, and write to file
    DISoutfile = open('data/{}.csv'.format(discordanceOUT), 'w')
    DISouthandle = ['accession,local_lineage,local_probability,local_pangoLEARN_version,local_status,seq_name,sample_col_date,sample_location,gisaid_lineage,gisaid_clade \n']
    for accession1, lineage1, prob, version, status, note, name, date, location, lineage2, clade in DIShandle:
        DISouthandle.append('{},{},{},{},{},{},{},{},{},{},{}\n'.format( accession1, lineage1, prob, version, status, note, name, date, location, lineage2, clade))
    DISoutfile.writelines(DISouthandle)

    DISoutfile.close()
    conn.close()


def parse_args():
    """ Command-line interface """
    parser = argparse.ArgumentParser(
        description="Update previous FASTA alignment with new sequences by"
                    "pairwise alignment to a reference genome."
    )
    parser.add_argument('--srcfile',
                        help='input, FASTA file with new sequences')
    parser.add_argument('--ref', default='data/NC_045512.fa',
                        help='Path to reference sequence.')
    parser.add_argument('--db', default='data/gsaid.db',
                        help='Name of database.')
    parser.add_argument('--sequencingmeta', '--sm',
                        help='tsv file from the sequencing meta table option on GSAID')
    parser.add_argument('--patientmeta', '--ptm',
                        help='tsv file from the patient meta data download option on GSAID.')
    parser.add_argument('--acknowledgementmeta', '--am',
                        help='XLS file containing acknowledgement data.')
    parser.add_argument('--outfasta', '-o',
                        help='Path to write outputfile for alignment')
    parser.add_argument('--targetdb', '-t',
                        help='Path to targetdb for migration')
    parser.add_argument('--lineagecsv', '-l',
                        help='Path to csv file containing pangolin output')
    parser.add_argument('--rawfasta',
                        help='Path to write outputfile for rawseqs')
    parser.add_argument('--lineageReport', action='store_true',
                        help='Export lineage (dis/con)cordance to LineageConcordance.csv, LineageDiscordance.csv')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    cursor, conn = open_connection(args.db)

    if args.srcfile is not None:
        iterate_fasta(args.srcfile)

    if args.sequencingmeta is not None:
        insert_sample_sequencing(cursor, args.sequencingmeta)

    if args.patientmeta is not None:
        insert_ptinfo(cursor, args.patientmeta)

    if args.acknowledgementmeta is not None:
        insert_acknowledgement(cursor, args.acknowledgementmeta)

    if args.outfasta is not None:
        export_fasta(cursor, args.outfasta)
        conn.close()

    if args.targetdb is not None:
        conn.close()
        migrate_entries(args.db, args.targetdb)

    if args.lineagecsv is not None:
        iterate_lineage_csv(cursor, args.lineagecsv)
        conn.commit()
        conn.close()
    if args.rawfasta is not None:
        dump_raw(args.rawfasta, db='data/gsaid.db')

    if args.lineageReport == True:
        dump_concordance(args.db, 'LineageConcordance', 'LineageDiscordance')
