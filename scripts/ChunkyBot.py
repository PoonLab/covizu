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
	old_dict = dict(convert_fasta(open(old_fasta)))
	new_dict = dict(convert_fasta(open(new_fasta)))
	val_diff={}
	miss_diff={}
	for key, value in old_dict.items():
		try:
			newval = new_dict[key]
			if newval != value:
				val_diff[key] = newval
		except:
			miss_diff[key] = value
	return miss_diff, val_diff



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
		'-o', '--outputdiffFasta', type=str, default='data/weeklydumps/diff.fasta',
		help='Directory for Output file containing sequences that are modified'
		)
	parser.add_argument(
		'-o2', '--outputmissFasta', type=str, default='data/weeklydumps/miss.fasta',
		help='Directory for Output file containing sequences that are missing'
		)
	return parser.parse_args()


if __name__ == '__main__':
	args = parse_args()
	download_folder = args.d.name
	driver = get_driver(download_folder=download_folder, executable_path=args.binpath)
	driver = login(driver=driver)
	bfiles = ['data/baseline/'+ s for s in os.listdir(args.baselinedir)]
	difffasta = {}
	missfasta = {}
	#predefined blocks one, two and three have <10k sequences
	srcfile = retrieve_genomes(driver=driver, start='2019-01-01', end='2020-04-16',
							   download_folder=download_folder)
	miss_dict, val_dict = compare_dicts(bfiles[0],srcfile)
	difffasta.update(val_dict)
	missfasta.update(miss_dict)
	time.sleep(300)

	srcfile = retrieve_genomes(driver=driver, start='2020-04-17', end='2020-05-08',
							   download_folder=download_folder)
	miss_dict, val_dict = compare_dicts(bfiles[1],srcfile)
	difffasta.update(val_dict)
	missfasta.update(miss_dict)
	time.sleep(300)

	srcfile = retrieve_genomes(driver=driver, start='2020-05-08', end='2020-05-15',
							   download_folder=download_folder)
	miss_dict, val_dict = compare_dicts(bfiles[2],srcfile)
	difffasta.update(val_dict)
	missfasta.update(miss_dict)
	time.sleep(300)

	#next do the chunks afterwards in weekly increments
	#find how many new chunks we need to get
	startdate = datetime(2020,5,16)
	today = datetime.now()
	delta = (today- startdate).days
	#find number of weekly blocks
	newchunks = delta//7 + 1


	for count in range(1, newchunks + 1):
		osdircount = 3
		#get the chunk dates
		enddate = startdate + timedelta(6)
		startmonthstr = str(startdate.month) if startdate.month > 9 else '0'+ str(startdate.month)
		startdaystr = str(startdate.day) if startdate.day > 9 else '0' + str(startdate.day)
		endmonthstr = str(enddate.month) if enddate.month > 9 else '0'+ str(enddate.month)
		enddaystr = str(enddate.day) if enddate.day > 9 else '0' + str(enddate.day)
		blockname = startmonthstr+startdaystr + '_' + endmonthstr+enddaystr

		srcfile = retrieve_genomes(driver=driver, start=datetime.strftime(startdate, '%Y-%m-%d'),
								   end=datetime.strftime(enddate, '%Y-%m-%d'), download_folder=download_folder)
		miss_dict, val_dict = compare_dicts(bfiles[osdircount],srcfile)
		difffasta.update(val_dict)
		missfasta.update(miss_dict)
		osdircount +=1
		time.sleep(300)
		startdate = enddate +timedelta(1)
		#set up for next iteration

	#move latest file directly to baseline
	shutil.move(srcfile, cwd+'/data/baseline/GISAID-'+blockname+ '.fasta')

	driver.quit()

	with open(args.outputdiffFasta, mode='w+') as out:
		for h,s in difffasta.items():
			out.write('>{}\n{}\n'.format(h, s))

	with open(args.outputmissFasta, mode='w+') as out2:
		for h,s in missfasta.items():
			out.write('>{}\n{}\n'.format(h, s))
