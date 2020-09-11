import os
import sys
import time
import re
import csv
from datetime import date, timedelta, datetime

from selenium import webdriver
from selenium.webdriver.firefox.options import Options

import shutil
import subprocess
import argparse
import tempfile
import getpass
import numpy as np

from gotoh2 import *
from autobot import get_driver, login, retrieve_genomes
from db_utils import pull_field, open_connection, insert_seq, find_seq, iterate_handle, insert_into_rawseqs, process_meta

todaystr = datetime.strftime(datetime.now(), '%Y-%m-%d')

def compare_fields(hash_dictionary, header_dictionary, new_fasta):
    """ Function to compare downloaded fasta to sequences present in database
    :params:
        :new_fasta: fasta file containing sequences downloaded from GISAID database
    :output:
        :modified seqs: list containing tuple [header, sequence]
    """
    new_fasta = convert_fasta(open(new_fasta))
    modified_seqs = {}
    for h,s in new_fasta:
        accession = h.split('|')[1]
        hashed_seq = hash(s.strip('N'))
        try:
            if hash_dictionary[accession] != hashed_seq or header_dictionary[accession] != h:
                #modified_seqs[h] = s
                pass
        except KeyError:
            modified_seqs[h] = s

    return modified_seqs

def retrieve_meta(driver, start, end, download_folder):
    """
    Retrieve meta data with a specified deposition date range in the GISAID database.
    Adding several time delays to avoid spamming the database.

    :param driver:  webdriver.Firefox object from login()
    :param start:  date in ISO format (yyyy-mm-dd)
    :param end:  date in ISO format (yyyy-mm-dd)
    :return:  paths to epi meta, patient meta
    """
    # find prefix variable
    element = driver.find_element_by_xpath("//div[@class='buttons container-slot']")
    htmlid_as_list = element.get_attribute('id').split('_')
    variable = htmlid_as_list[1]

    # trigger selection change
    time_string = '[id^="ce_' + variable + '"][id$="_input"]'

    driver.execute_script("document.querySelectorAll('{}')[2].value = '{}'".format(time_string, start))
    driver.execute_script("document.querySelectorAll('{}')[2].onchange()".format(time_string))

    driver.execute_script("document.querySelectorAll('{}')[3].value = '{}'".format(time_string, end))
    driver.execute_script("document.querySelectorAll('{}')[3].onchange()".format(time_string))

    driver.execute_script("document.querySelectorAll('[id^=\"ce_{}\"][id$=_input]')[2].onchange()".format(variable))
    time.sleep(15)

    print('selecting all seqs')
    element = driver.find_element_by_xpath("//*[contains(text(), 'Total')]")
    count = element.get_attribute('innerHTML').split()[1].replace(',', '')
    if int(count) > 10000:
        time.sleep(15)

    checkbox = driver.find_element_by_xpath("//span[@class='yui-dt-label']/input[@type='checkbox']")
    checkbox.click()
    time.sleep(5)

    #download two meta data files - epi, tech
    metafiles = []
    metafiles.append(wait_for_download(driver, download_folder, 'meta_epi'))
    metafiles.append(wait_for_download(driver, download_folder, 'meta_tech'))

    #reset selection
    element = driver.find_element_by_xpath("//button[@class='sys-event-hook sys-form-button']")
    element.click()

    return metafiles

def wait_for_download(driver, download_folder, target):
    # download pt status
    element = driver.find_element_by_xpath("//*[contains(text(), 'Download')]")
    driver.execute_script("arguments[0].click();", element)
    time.sleep(5)
    # switch to iframe to download
    driver.switch_to.frame(driver.find_element_by_tag_name("iframe"))
    print("Downloading {}.".format(target))
    time.sleep(5)
    file_query = "input[type='radio'][value='{}']".format(target)
    radio = driver.find_element_by_css_selector(file_query).click()

    button = driver.find_element_by_xpath(
        "//*[contains(text(), 'Download')]//ancestor::div[@style='float: right']"
    )
    button.click()
    time.sleep(5)

    # wait for download to complete
    while True:
        files = os.listdir(download_folder)
        if len(files) == 0:
            # download has not started yet
            time.sleep(10)
            continue
        if any([f.endswith('.part') for f in files]):
            time.sleep(5)
            continue
        break

    # reset browser
    time.sleep(60)
    driver.switch_to.default_content()
    #Get newest file in directory
    downloaded_dir= os.listdir(download_folder)
    abs_dir= []
    for path in downloaded_dir:
        abs_dir.append(os.path.join(download_folder, path))
    downloaded_file = max(abs_dir, key=os.path.getctime)
    print('Download of {} complete.'.format(downloaded_file))

    return downloaded_file

def parse_args():
    """ Command line interface """
    parser = argparse.ArgumentParser(
        description="Python3 Script to Automate retrieval of genomes deposited in a given day."
    )
    parser.add_argument(
        '-d', '--dir', type=str, default=tempfile.TemporaryDirectory(),
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
    parser.add_argument(
        '--debug', action='store_true', help='Saves missing sequences into missing.fa flat file'
    )
    return parser.parse_args()


if __name__ == '__main__':
    #initialize webdriver
    args = parse_args()
    download_folder = args.dir.name
    driver = get_driver(download_folder=download_folder, executable_path=args.binpath)
    driver = login(driver=driver)
    log = 'Changed headers \n' #debug tool

    _, refseq = convert_fasta(open(args.ref))[0]

    #load in existing file names with their pre-defined date ranges
    bdates= [['2019-01-01','2020-04-16'],['2020-04-17','2020-05-08'],['2020-05-09','2020-05-15'], ['2020-05-16', '2020-05-31'], ['2020-06-01','2020-06-10'], \
    ['2020-06-11','2020-06-22'], ['2020-06-23','2020-07-04'], ['2020-07-05','2020-07-21'], ['2020-07-22','2020-08-07'], ['2020-08-08','2020-08-24'], \
    ['2020-08-25','2020-09-05']]

    #find the date ranges for all non-predefined chunks
    startdate = datetime(2020, 9, 6)
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

    #initialize dictionaries to compare with
    cursor, conn = open_connection(args.database)
    hash_dictionary = pull_field(cursor, 'unaligned_hash')
    header_dictionary= pull_field(cursor, 'header')
    conn.close()
    modified_seqs= {}

    missing_from_gsaid = list(header_dictionary.keys())

    #download all chunks, and compare sequences within
    for start, end in bdates:
        srcfile = retrieve_genomes(driver=driver, start=start, end=end, download_folder=download_folder)
        print(srcfile)
        #update modified_seqs dictionary
        modified_seqs.update(compare_fields(hash_dictionary, header_dictionary, srcfile))
        raw_accessions = insert_into_rawseqs(args.database, srcfile)
        #retrieve meta data
        time.sleep(60)
        metafiles = retrieve_meta(driver, start=start, end=end, download_folder=download_folder)
        process_meta(args.database, metafiles)
        time.sleep(60)
        missing_from_gsaid= np.setdiff1d(missing_from_gsaid, raw_accessions)
    if args.debug:
        #:DEBUG:
        debugout= open('missing.fa', 'w')
        for h,s in modified_seqs.items():
            debugout.write('>{}\n{}\n'.format(h,s))
        debugout.close()

    driver.quit()
    print('Number of seqs removed from GISAID database: {}\n'.format(len(missing_from_gsaid)))
    #call the updater, passing modified sequences as a list
    iterate_handle(modified_seqs.items(), args.ref, database = args.database)
