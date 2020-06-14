import os
import sys
import time
from datetime import date, timedelta, datetime

from selenium import webdriver
from selenium.webdriver.firefox.options import Options

import shutil
import subprocess
import argparse
import tempfile
import getpass

from gotoh2 import *
from autobot import get_driver, login, retrieve_genomes


todaystr = datetime.strftime(datetime.now(), '%Y-%m-%d')
cwd = os.getcwd()

def compare_dicts(old_fasta, new_fasta):
	""" Function to compare contents of two fasta files
	:params:
		:old_fasta: fasta file containing sequences already downloaded
		:new_fasta: fasta file containing sequences downloaded from GISAID database
	:output:
		:new_old_header: list containing new header and existing header in local files
		:seq_diff: dictionary containing header:sequence pairs that need to be replaced
	"""
	new_old_header = []
	seq_diff = {}

	old_dict = dict(convert_fasta(open(old_fasta)))
	new_dict = dict(convert_fasta(open(new_fasta)))

	#create dictionaries with key:val pair being ascension no : header
	old_header_dict = {}
	for header in old_dict.keys():
		old_header_dict[header.split('|')[1]] = header
	new_headers_dict = {}
	for header in new_dict.keys():
		new_headers_dict[header.split('|')[1]] = header

	#compare header dicts
	for a, h in old_header_dict.items():
		try:
			newheader = new_headers_dict[a]
			if newheader != h:
				new_old_header.append((newheader, h))
		except:
			pass
	#compare sequences
	for h, s in old_dict.items():
		try:
			newseq = new_dict[h]
			if newseq != s:
				seq_diff[h] = newval
		except:
			pass
	return new_old_header, seq_diff

def download_compare(baselinefile, startdate, endate, download_folder=download_folder, driver=driver):
	""" Wrapper function for download & compare function
	:params:
		:baselinefile: string, path to baseline file
		:startdate: string, iso date (yyyy-mm-dd) for start of chunk
		:enddate: string, iso date (yyyy-mm-dd) for end of chunk
		\\optional\
		:download_folder: tempfolder where seqs are to be downloaded to
		:driver: webdriver

	:output:
		:miss_diff: dictionary containing key
		:new_old_header: list containing tuples of headers that need to be modified

	"""
	srcfile = retrieve_genomes(driver=driver, start=startdate, end=endate,
                                                           download_folder=download_folder)
	#update baseline files if this is the last chunk
	if baselinefile == bfiles[-1]:
		shutil.move(srcfile, cwd+ '/' +baselinefile)

	new_old_header, seq_diff= compare_dicts(baselinefile,srcfile) #compare & report needed changes
	modify_fasta(baselinefile, new_old_header, seq_diff) #change the baselinefile with changes

	return new_old_header, seq_diff

def modify_fasta(fasta, new_old_header, seq_diff):
	"""Function that replaces headers and sequences
	:params:
		:fasta: str, path to fasta file to be modified
		:new_old_header: list containing tuples of headers that need to be modified
		:seq_diff: dictionary containing sequences that need to be changed
	:output:
	"""
	pass
def parse_args():
	""" Command line interface """
	parser = argparse.ArgumentParser(
		description="Python3 Script to Automate retrieval of genomes deposited in a given day."
	)
	parser.add_argument(
		'-d', '-dir', type=str, default=tempfile.TemporaryDirectory(),
		help="Temporary directory to download files."
	)
	parser.add_argument(
		'-b', '--binpath', type=str, default='/usr/local/bin/geckodriver',
		help='Path to geckodriver binary executable'
	)
	parser.add_argument(
		'-l', '--baselinedir', type=str, default='data/baseline',
		help='Folder containing baseline downloads'
	)
	parser.add_argument(
		'-a', '--alignment', type=str, default='data/gisaid-aligned.fa',
		help='Path to alignment in Fasta'
	)
	return parser.parse_args()


if __name__ == '__main__':
	#initialize webdriver
	args = parse_args()
	download_folder = args.d.name
	driver = get_driver(download_folder=download_folder, executable_path=args.binpath)
	driver = login(driver=driver)
	log = 'Changed headers \n' #debug tool


	#load in existing file names with their pre-defined date ranges
	bfiles = ['data/baseline/'+ s for s in os.listdir(args.baselinedir)]
	bdates= [('2019-01-01','2020-04-16'),('2020-04-17','2020-05-08'),('2020-05-08','2020-05-15')]
	difffasta = {}
	missfasta = {}

	#find the date ranges for all non-predefined chunks
	startdate = datetime(2020, 5, 16)
	today = datetime.now()
	delta = (today- startdate).days
	newchunks = delta//7 + 1 	#find number of weekly blocks
	#Calculate the date ranges for each chunk
	for count in range(1, newchunks + 1):
		enddate = startdate + timedelta(6)
		startmonthstr = str(startdate.month) if startdate.month > 9 else '0'+ str(startdate.month)
		startdaystr = str(startdate.day) if startdate.day > 9 else '0' + str(startdate.day)
		endmonthstr = str(enddate.month) if enddate.month > 9 else '0'+ str(enddate.month)
		enddaystr = str(enddate.day) if enddate.day > 9 else '0' + str(enddate.day)
		blockname = 'data/baseline/'+ startmonthstr+startdaystr + '_' + endmonthstr+enddaystr
			bfiles.append(blockname)
		bdates.append(( datetime.strftime(startdate, '%Y-%m-%d'),
				datetime.strftime(enddate, '%Y-%m-%d')
			     ))
		startdate = enddate +timedelta(1)

	#start downloading & comparing
	for file, index in enumerate(bfiles):
		new_old_header, seq_diff = download_compare(file, bdates[index][0], bdates[index][1],
							    download_folder=download_folder, driver=driver)
		modify_alignment(new_old_header, seq_diff, gsaid_alignment)
		for tuple in new_older_header:
			log+= tuple + '\n'
		for header in seq_diff.keys():
			log+= header +'\n'


	driver.quit()

	#debug section
	print(log)
