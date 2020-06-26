import os
import sys
import time
import re
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
from db_utils import pull_field, open_connection, insert_seq, find_seq

todaystr = datetime.strftime(datetime.now(), '%Y-%m-%d')

def compare_fields(hash_dictionary, header_dictionary, new_fasta):
    """ Function to compare downloaded fasta to sequences present in database
    :params:
        :new_fasta: fasta file containing sequences downloaded from GISAID database
    :output:
        :modified seqs: list containing tuple [header, sequence]
    """
    new_fasta = convert_fasta(open(new_fasta))
    modified_seqs = []

    for h,s in new_fasta:
        accession = h.split('|')[1]
        hashed_seq = hash(s.strip('N'))
        try:
            if hash_dictionary[accession] != hashed_seq or header_dictionary[accession] != h:
                modified_seqs.append([h,s])
        except KeyError:
            modified_seqs.append([h,s])

    return modified_seqs

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
    parser.add_argument(
        '-db', '--database', type=str, default = 'data/gsaid.db',
        help='Path to sqlite3 database containing sequences'
    )
    parser.add_argument(
        '-r', '--ref', type=str, default = 'data/NC_045512.fa',
        help='Path to ref seq'
    )
    return parser.parse_args()


if __name__ == '__main__':
    #initialize webdriver
    args = parse_args()
    download_folder = args.d.name
    driver = get_driver(download_folder=download_folder, executable_path=args.binpath)
    driver = login(driver=driver)
    log = 'Changed headers \n' #debug tool

    _, refseq = convert_fasta(open(args.ref))[0]

    #load in existing file names with their pre-defined date ranges
    bdates= [['2019-01-01','2020-04-16'],['2020-04-17','2020-05-08'],['2020-05-08','2020-05-15']]

    #find the date ranges for all non-predefined chunks
    startdate = datetime(2020, 5, 16)
    today = datetime.now()
    delta = (today- startdate).days
    newchunks = delta//7 + 1    #find number of weekly blocks
    #Calculate the date ranges for each chunk
    for count in range(1, newchunks):
        enddate = startdate + timedelta(6)
        startmonthstr = str(startdate.month) if startdate.month > 9 else '0'+ str(startdate.month)
        startdaystr = str(startdate.day) if startdate.day > 9 else '0' + str(startdate.day)
        endmonthstr = str(enddate.month) if enddate.month > 9 else '0'+ str(enddate.month)
        enddaystr = str(enddate.day) if enddate.day > 9 else '0' + str(enddate.day)
        bdates.append([ datetime.strftime(startdate, '%Y-%m-%d'),
                datetime.strftime(enddate, '%Y-%m-%d')]
                 )
        startdate = enddate +timedelta(1)

    #start downloading & comparing

    cursor, conn = open_connection(args.database)

    hash_dictionary = pull_field(cursor, 'unaligned_hash')
    header_dictionary= pull_field(cursor, 'header')
    modified_seqs= []

    for start, end in bdates:
        srcfile = retrieve_genomes(driver=driver, start=start, end=end, download_folder=download_folder)
        modified_seqs+= compare_fields(hash_dictionary, header_dictionary, srcfile)

    driver.quit()
    print(len(modified_seqs))
    print(modified_seqs)

    for h,s in modified_seqs:
        aligned = find_seq(conn, s, refseq)
        insert_seq(cursor, s, h, aligned)

    conn.commit()
    conn.close()


    #debug section
    #:TODO
