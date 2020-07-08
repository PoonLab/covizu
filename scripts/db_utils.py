import sqlite3
import gotoh2
import argparse
import threading
import queue
import time


class AlignerThread(threading.Thread):
    def __init__(self, in_queue, out_queue):
        threading.Thread.__init__(self)
        self.in_queue = in_queue
        self.out_queue = out_queue

    def run(self):
        # unpack data payload
        header, seq, database, refseq = self.in_queue.get()
        # open db to check if sequences already aligned
        cursor, conn = open_connection(database)
        alignedseq = find_seq(conn, seq, refseq)
        if alignedseq == 0:
            self.in_queue.task_done()
            return 1
        self.out_queue.put((header, seq, alignedseq))
        conn.close()
        self.in_queue.task_done()
        return 0


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
                 'PRIMARY KEY, header VARCHAR(255), unaligned_hash VARCHAR(100), ' \
                 'aligned BLOB);'
    cur.execute(seqs_table)

    sample_table = "CREATE TABLE IF NOT EXISTS SAMPLE (accession VARCHAR(255) " \
                   "PRIMARY KEY, virus_name VARCHAR(255), collection_date VARCHAR(10), " \
                   "location VARCHAR(30));"
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

    conn.commit()
    return cur, conn


def insert_seq(cursor, seq, header, aligned):
    """ Wrapper function for aligning & inserting into SEQUENCES table
        :params:
            :cursor: sqlite database handler
            :seq: raw sequence
            :header: fasta header
    """
    vars= [header.split('|')[1], header, hash(seq.strip('N')), aligned]
    result = cursor.execute("REPLACE INTO SEQUENCES('accession', 'header', 'unaligned_hash', 'aligned') VALUES(?, ?, ?, ?)", vars)
    return 0


def find_seq(conn, seq, refseq):
    """
    Function that either returns aligned seq from database or performs procrust align itself
    :params:
        :cursor: sqlite db handler
        :seq: raw sequence
        :ref: NC_045512 sequence
    """
    if len(seq) < 5000:
        return 0

    # hash the raw sequence to find the hashkey in aligned table
    rawhash = hash(seq.strip('N'))
    result = conn.cursor().execute('SELECT aligned_seq FROM HASHEDSEQS WHERE `hash_key` = ?',
                                   [rawhash]).fetchall()
    if len(result) == 1:
        aligned = result[0][0]
    else:
        aligner = gotoh2.Aligner()
        aligned = gotoh2.procrust_align(refseq, seq, aligner)[0]
        vars = [hash(seq.strip('N')), aligned]
        conn.cursor().execute('REPLACE INTO HASHEDSEQS(`hash_key`, '
                              '`aligned_seq`) VALUES(?,?);', vars)
        conn.commit()

    return aligned


def insert_sample_sequencing(cursor, tsvFile):
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
        handle = handle.readlines()[2:]
        for row in handle.split('	'):
            sample_vars = [row[1], row[0], row[2], row[3]]
            result1 = cursor.execute("REPLACE INTO SAMPLE('accession', 'virus_name', "
                                     "'collection_date', 'location') VALUES(?, ?, ?, ?)", sample_vars)
            sequencing_vars = [row[1], row[8], row[9], row[6]]
            result2 = cursor.execute("REPLACE INTO SEQUENCING('accession', "
                                     "'platform', 'assembly_method', 'sample_type') "
                                     "VALUES(?, ?, ?, ?)", sequencing_vars)
            conn.commit()


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
        handle = handle.readlines()[3:]
        for row in handle.split('	'):
            vars = [row[0], row[4], row[5], row[6]]
            result = cursor.execute("REPLACE INTO ACKNOWLEDGEMENTS('accession', 'originating_lab', 'sub_lab', "
                                    "'authors') VALUES(?, ?, ?, ?)", vars)
            conn.commit()


def insert_ptinfo(cursor, tsvFile):
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
        handle = handle.readlines()[2:]
        for row in handle.split('	'):
            vars = [row[1], row[5], row[6], row[7]]
            result = cursor.execute("REPLACE INTO PTINFO('accession', 'gender', 'age', "
                                    "'ptstatus') VALUES(?, ?, ?, ?)", vars)
            conn.commit()


def pull_field(cursor, field):
    """
    Wrapper function to return dictionary key-val pair being accession-field
    """
    field = ["`" + field + "`"]
    result = cursor.execute("SELECT `accession`, ? FROM SEQUENCES;", field)
    field_dict= dict(result.fetchall())
    return field_dict


def iterate_fasta_threaded(fasta, ref, database = 'data/gsaid.db'):
    """
    Function to iterate through fasta file
    :param cursor: sqlite database handler
    :param fasta: file containing sequences
    :param ref: file containing reference sequence, default-> NC_04552.fa
    """
    handle = gotoh2.iter_fasta(open(fasta))
    _, refseq = gotoh2.convert_fasta(open(ref))[0]

    in_queue = queue.Queue()
    out_queue = queue.Queue()

    # stream records into queue
    for h, s in handle:
        in_queue.put((h, s, database, refseq))

    # pairwise align sequences in queue
    for seq in range(0,in_queue.qsize()):
        AlignerThread(in_queue, out_queue).start()

    # block until queue is processed
    in_queue.join()

    # update database with aligned sequences
    cursor, conn = open_connection(database)
    while not out_queue.empty():
        header, seq, alignedseq = out_queue.get()
        insert_seq(cursor, seq, header, alignedseq)
        conn.commit()
        out_queue.task_done()
    conn.close()
    out_queue.join()

def iterate_fasta(fasta, ref, database = 'data/gsaid.db'):
    """
    Function to iterate through fasta file *non threaded*
    :param cursor: sqlite database handler
    :param fasta: file containing sequences
    :param ref: file containing reference sequence, default-> NC_04552.fa
    """
    handle = gotoh2.convert_fasta(open(fasta))
    _, refseq = gotoh2.convert_fasta(open(ref))[0]
    cursor, conn = open_connection(database)


    # stream records into queue
    for header, seq in handle:
        alignedseq = find_seq(conn, seq, refseq)
        insert_seq(cursor, seq, header, alignedseq)
        conn.commit()
    conn.close()


def iterate_handle(handle, ref, database = 'data/gsaid.db'):
    """
    Function to iterate through list [(header,seq)....] passed by ChunkyBot ##NON THREADED##
    :param cursor: sqlite database handler
    :param handle: list containing tuples (header, raw seq)
    :param ref: file containing reference sequence, default-> NC_04552.fa
    """

    _, refseq = gotoh2.convert_fasta(open(ref))[0]
    cursor, conn = open_connection(database)

    #align and insert into database
    for header, seq in handle:
        alignedseq = find_seq(conn, seq, refseq)
        #if alignedseq returns 0; raw seq length <5,000
        if alignedseq == 0:
            print('Sequence {} with a length of {} cannot be aligned.'.format(header, len(seq)))
        else:
            print('Aligning {}.'.format(header))
            insert_seq(cursor, seq, header, alignedseq)
            conn.commit()

def iterate_handle_threaded(handle, ref, database = 'data/gsaid.db'):
    """
    Function to iterate through handle passed by Chunky bot
    :param cursor: sqlite database handler
    :param handle: list containing tuples (header, raw seq)
    :param ref: file containing reference sequence, default-> NC_04552.fa
    """

    _, refseq = gotoh2.convert_fasta(open(ref))[0]

    in_queue = queue.Queue()
    out_queue = queue.Queue()

    # stream records into queue
    for h, s in handle:
        in_queue.put((h, s, database, refseq))

    # pairwise align sequences in queue
    for seq in range(0,in_queue.qsize()):
        AlignerThread(in_queue, out_queue).start()

    # block until queue is processed
    in_queue.join()

    # update database with aligned sequences
    cursor, conn = open_connection(database)
    while not out_queue.empty():
        header, seq, alignedseq = out_queue.get()
        insert_seq(cursor, seq, header, alignedseq)
        conn.commit()
        out_queue.task_done()
    conn.close()
    out_queue.join()

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
    #get existing accessions in target db, create list to derive missing accesions
    existing_accessions = list(oldcursor.execute('SELECT accession from SEQUENCES;').fetchall())
    EA_list = []
    for accession in existing_accessions:
        EA_list.append(accession[0])

    #get all sequences from the source db and insert them into target
    All_chunk_seqs = list(chunkcursor.execute('SELECT * FROM SEQUENCES;').fetchall())
    count = 0
    for accession, header, unaligned_hash, aligned in All_chunk_seqs:
        if accession not in EA_list:
            count+=1
            oldcursor.execute('INSERT INTO SEQUENCES (accession, header, unaligned_hash, aligned) VALUES (?,?,?,?);', [accession, header, unaligned_hash, aligned])
            oldconnect.commit()

    print('Updated {}, inserted {} seqs.'.format(old_db, str(count)))

    #get all existing hashes from target db, create list to derive missing accesions
    existing_hashes = list(oldcursor.execute('SELECT hash_key FROM HASHEDSEQS').fetchall())
    EH_list = []
    for hash_key in existing_hashes:
        EH_list.append(hash_key[0])

    hash_count = 0
    All_chunk_hashes = list(chunkcursor.execute('SELECT * FROM HASHEDSEQS;').fetchall())
    for hash_key, aligned_seq in All_chunk_hashes:
        if hash_key not in EH_list:
            hash_count+=1
            oldcursor.execute('INSERT INTO HASHEDSEQS (hash_key, aligned_seq) VALUES (?,?)', [hash_key, aligned_seq])
            oldconnect.commit()
    print('Updated {}, inserted {} hashes.'.format(old_db, str(hash_count)))

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
    parser.add_argument('--targetdb',
                        help='Name of target database')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    conn, cursor = open_connection(args.db)

    if args.srcfile is not None:
        iterate_fasta(args.srcfile, args.ref)

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


