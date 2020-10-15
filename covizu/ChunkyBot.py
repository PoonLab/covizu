import os
import sys
import time
import re
import csv
import datetime

from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.common.exceptions import NoSuchElementException

import shutil
import subprocess
import argparse
import tempfile
import getpass
import numpy as np

from covizu.autobot import *
from covizu.pangorider import *
from covizu.db_utils import *

todaystr = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d')

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
    if metafiles[0] == 0:
        return []
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
    time.sleep(10)

    file_query = "input[type='radio'][value='{}']".format(target)

    try:
        radio = driver.find_element_by_css_selector(file_query).click()
    except:
        return 0 

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

def check_for_changes(srcfile, database):
    """
    Iterates through srcfile, compares headers and raw sequences
    Returns sequence tuple list that contains sequences that need to be re-aligned
    """

    fasta = convert_fasta(open(srcfile, 'r'))
    sequence_handle = []

    for h,s in fasta:
        result = report_changes(database, h,s)
        if result == 1:
            print('New header inserted for {}'.format(h))
        if result == 0:
            continue
        else:
            sequence_handle.append(result)

    return sequence_handle




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
        '-r', '--ref', type=str, default = 'data/LR757995.fa',
        help='Path to ref seq'
    )
    parser.add_argument(
        '--debug', action='store_true', help='Saves missing sequences into missing.fa flat file'
    )
    parser.add_argument('--filterout', default = 'debug/filtered.log', type =argparse.FileType('w+'),
        help='Log for filtered seqs')

    parser.add_argument('--indicies', help='Indicies to keep, defaults to [265:29674]')

    site_packages = next(p for p in sys.path if 'site-packages' in p)
    parser.add_argument('--pangolindir', help ='PangoLEARN data dir, defaults to .../site_packages/panoLEARN/data/',
        default=site_packages+'/pangoLEARN/data/')
    return parser.parse_args()


if __name__ == '__main__':
    #initialize webdriver
    args = parse_args()
    download_folder = args.dir.name
    driver = get_driver(download_folder=download_folder, executable_path=args.binpath)
    driver = login(driver=driver)
    log = 'Changed headers \n' #debug tool

    _, refseq = convert_fasta(open(args.ref))[0]

    #load in existing file names with their pre-defined date ranges, #2020-09-28 has 11k sequences
    #removing blocks before July ['2019-01-01', '2020-04-16'], ['2020-04-17', '2020-05-08'], ['2020-05-09', '2020-05-15'], ['2020-05-16', '2020-05-31'], ['2020-06-01', '2020-06-10'], ['2020-06-11', '2020-06-22'],

    bdates= [ ['2020-06-23', '2020-07-04'], ['2020-07-05', '2020-07-21'], ['2020-07-22', '2020-08-07'], ['2020-08-08', '2020-08-24'], \
        ['2020-08-25', '2020-09-05'], ['2020-09-06', '2020-09-12'], ['2020-09-13', '2020-09-19'], ['2020-09-20', '2020-09-27'], ['2020-09-29', '2020-10-03']]


    #find the date ranges for all non-predefined chunks
    startdate = datetime.datetime(2020, 10, 4)
    today = datetime.datetime.now()
    delta = (today- startdate).days
    newchunks = delta//7 + 1    #find number of weekly blocks
    #Calculate the date ranges for each chunk
    for count in range(1, newchunks):
        enddate = startdate + timedelta(6)
        startmonthstr = str(startdate.month) if startdate.month > 9 else '0'+ str(startdate.month)
        startdaystr = str(startdate.day) if startdate.day > 9 else '0' + str(startdate.day)
        endmonthstr = str(enddate.month) if enddate.month > 9 else '0'+ str(enddate.month)
        enddaystr = str(enddate.day) if enddate.day > 9 else '0' + str(enddate.day)
        bdates.append([ datetime.datetime.strftime(startdate, '%Y-%m-%d'),
                datetime.datetime.strftime(enddate, '%Y-%m-%d')]
                 )
        startdate = enddate +timedelta(1)

    #start downloading & comparing

    #initialize dictionaries to compare with
    cursor, conn = open_connection(args.database)
    header_dictionary= pull_field(cursor, 'header')
    conn.close()
    modified_seqs= {}
    sequence_handle = []

    missing_from_gsaid = list(header_dictionary.keys())

    #download all chunks, and compare sequences within
    chunk_count = 0 # reset browser after 10 chunks

    for start, end in bdates:
        srcfile = retrieve_genomes(driver=driver, start=start, end=end, download_folder=download_folder)
        print(srcfile)
        raw_accessions = insert_into_rawseqs(args.database, srcfile)
        sequence_handle += check_for_changes(srcfile, args.database)
        #if sequence_handle gets too large, run it through the classifer to reduce memory usage
        if len(sequence_handle) > 800:
            filtered_handle = filter_seqs(sequence_handle, args.filterout, max_prop_n=0.05, minlen=29000)
            classify_and_insert(args.pangolindir+ 'decisionTreeHeaders_v1.joblib', args.pangolindir+'decisionTree_v1.joblib', filtered_handle, args.indicies, args.database)
            sequence_handle = []
        #retrieve meta data
        time.sleep(60)
        metafiles = retrieve_meta(driver, start=start, end=end, download_folder=download_folder)
        print(start, end)
        if len(metafiles) == 0:
            continue
        process_meta(args.database, metafiles)
        time.sleep(60)
        missing_from_gsaid= np.setdiff1d(missing_from_gsaid, raw_accessions)
        #reset browser & browser counter after hitting 10 chunks
        if chunk_count == 10:
            driver.quit()
            driver = get_driver(download_folder=download_folder, executable_path=args.binpath)
            driver = login(driver=driver)
            chunk_count = 0
        chunk_count+=1

    if args.debug:
        #:DEBUG:
        debugout= open('missing.fa', 'w')
        for h,s in modified_seqs.items():
            debugout.write('>{}\n{}\n'.format(h,s))
        debugout.close()

    driver.quit()
    print('Number of seqs removed from GISAID database: {}\n'.format(len(missing_from_gsaid)))

    #handle the remaining sequences that are missing
    filtered_handle = filter_seqs(sequence_handle, args.filterout, max_prop_n=0.05, minlen=29000)
    if len(filtered_handle) > 0:
        classify_and_insert(args.pangolindir+ 'decisionTreeHeaders_v1.joblib', args.pangolindir+'decisionTree_v1.joblib', filtered_handle, args.indicies, args.database)
    else:
        print('No sequences to update')
