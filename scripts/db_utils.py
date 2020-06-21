import sqlite3
import gotoh2
import argparse
import threading

#db utils

class myThread(threading.Thread):

	def __init__(self, header, seq, database, refseq):
		threading.Thread.__init__(self)
		self.header = header
		self.seq = seq
		self.refseq = refseq
		self.database = database

	def run(self):
		#open db to check if sequences already aligned
		self.cursor, self.conn = open_connection(self.database)
		self.aligned = find_seq(self.conn, self.seq, self.refseq)
		#reopen db to insert sequence
		self.cursor, self.conn = open_connection(self.database)
		insert_seq(self.cursor, self.seq, self.header, self.alignedseq)
		self.conn.close()


def open_connection(database):
	""" open connection to database, initialize tables if they don't exist
		:params:
			:database: str, name of database
		:out:
			:cursor: interactive sql object containing tables
	"""
	conn = sqlite3.connect(database, check_same_thread=False)
	cur = conn.cursor()

	#create tables if they don't exist
	seqs_table = 'CREATE TABLE IF NOT EXISTS SEQUENCES (accession VARCHAR(255) PRIMARY KEY, header VARCHAR(255), unaligned BLOB, aligned BLOB);' #sequences
	cur.execute(seqs_table)
	sample_table= "CREATE TABLE IF NOT EXISTS SAMPLE (accession VARCHAR(255) PRIMARY KEY, virus_name VARCHAR(255), collection_date VARCHAR(10), location VARCHAR(30));"
	cur.execute(sample_table)
	sequencing_table =  "CREATE TABLE IF NOT EXISTS SEQUENCING (accession VARCHAR(255) PRIMARY KEY, platform VARCHAR(100), assembly_method VARCHAR(100), sample_type VARCHAR(100))"
	cur.execute(sequencing_table)
	acknowledgement_table = "CREATE TABLE IF NOT EXISTS ACKNOWLEDGEMENTS (accession VARCHAR(255) PRIMARY KEY, originating_lab VARCHAR(100), sub_lab VARCHAR(100), authors BLOB)"
	cur.execute(acknowledgement_table)
	ptinfo_table = "CREATE TABLE IF NOT EXISTS PTINFO (accession VARCHAR(255) PRIMARY KEY, gender VARCHAR(30), age INT, ptstatus VARCHAR(20))"
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
	vars= [header.split('|')[1], header, seq, aligned]

	result = cursor.execute("REPLACE INTO SEQUENCES('accession', 'header', 'unaligned', 'aligned') VALUES(?, ?, ?, ?)", vars)
	conn.commit()


def find_seq(conn, seq, refseq):
	""" Function that either returns aligned seq from database or performs procrust align itself
		:params:
			:cursor: sqlite db handler
			:seq: raw sequence
			:ref: NC_045512 sequence
	"""
	rawhash = hash(seq) #hash the raw sequence to find the hashkey in aligned table
	result = conn.cursor().execute('SELECT aligned_seq FROM HASHEDSEQS WHERE `hash_key` = ?', [rawhash]).fetchall()
	conn.close() #either way close the connection after running first query

	if len(result) == 1:
		aligned =  result[0][0]
	else:
		aligner = gotoh2.Aligner()
		aligned = gotoh2.procrust_align(refseq, seq, aligner)[0]

	return aligned

def insert_sample_sequencing(cursor, tsvFile):
	""" Wrapper function for inserting into SAMPLE & SEQUENCING table
		#FROM SEQUENCING TECH METADATA
		:params:
			:cursor: sqlite database handler
			:tsvFile: tabs separated values file containing data, from the  Sequencing technology metadata download option
		:out:
			:debug: report
	"""

	with open(tsvFile) as handle:
		handle = handle.readlines()[2:]
		for row in handle.split('	'):
			sample_vars = [row[1], row[0], row[2], row[3]]
			result1 = cursor.execute("REPLACE INTO SAMPLE('accession', 'virus_name', 'collection_date', 'location') VALUES(?, ?, ?, ?)", sample_vars)
			sequencing_vars = [row[1], row[8], row[9], row[6]]
			result2 = cursor.execute("REPLACE INTO SEQUENCING('accession', 'platform', 'assembly_method', 'sample_type') VALUES(?, ?, ?, ?)", sequencing_vars)
			conn.commit()

def insert_acknowledgement(cursor,tsvFile):
	""" Wrapper function for processing tsvFile and inserting into PTINFO table
		#FROM Patient status metadata
		:params:
			:cursor: sqlite database handler
			:tsvFile: tabs separated values file containing data, from the  Sequencing technology metadata download option
		:out:
			:debug: report
	"""


	with open(tsvFile) as handle:
		handle = handle.readlines()[3:]
		for row in handle.split('	'):
			vars = [row[0], row[4], row[5], row[6]]
			result = cursor.execute("REPLACE INTO ACKNOWLEDGEMENTS('accession', 'originating_lab', 'sub_lab', 'authors') VALUES(?, ?, ?, ?)", vars)
			conn.commit()

def insert_ptinfo(cursor, tsvFile):
	""" Wrapper function for processing tsvFile and inserting into PTINFO table
		#FROM Patient status metadata
		:params:
			:cursor: sqlite database handler
			:tsvFile: tabs separated values file containing data, from the patient status metadata download option
		:out:
			:debug: report
	"""

	with open(tsvFile) as handle:
		handle = handle.readlines()[2:]
		for row in handle.split('	'):
			vars = [row[1], row[5], row[6], row[7]]
			result = cursor.execute("REPLACE INTO PTINFO('accession', 'gender', 'age', 'ptstatus') VALUES(?, ?, ?, ?)", vars)
			conn.commit()

def iterate_fasta(fasta, ref, database = 'data/gsaid.db'):
	""" Function to iterate through fasta file
		:params:
			:cursor: sqlite database handler
			:fasta: file containing sequences
			:ref: file containing reference sequence, default-> NC_04552.fa
	"""
	handle = gotoh2.iter_fasta(open(fasta))
	_, refseq = gotoh2.convert_fasta(open(ref))[0]
	for h, s in handle:
		myThread(h, s, database, refseq).start()
		#insert_seq(cursor, s, h, refseq)


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
	parser.add_argument('--db', default = 'data/gsaid.db',
			help='Name of database.')
	parser.add_argument('--sequencingmeta', '--sm',
			help='tsv file from the sequencing meta table option on GSAID')
	parser.add_argument('--patientmeta', '--ptm',
			help='tsv file from the patient meta data download option on GSAID.')
	parser.add_argument('--acknowledgementmeta', '--am',
			help='XLS file containing acknowledgement data.')
	return parser.parse_args()

if __name__ == '__main__':
	args= parse_args()

	if args.srcfile is not None:
		iterate_fasta(cursor, args.srcfile, args.ref)

	if args.sequencingmeta is not None:
		insert_sample_sequencing(cursor, args.sequencingmeta)

	if args.patientmeta is not None:
		insert_ptinfo(cursor, args.patientmeta)

	if args.acknowledgementmeta is not None:
		insert_acknowledgement(cursor, args.acknowledgementmeta)





