import os
import sys
import time
from datetime import date, timedelta

from selenium import webdriver
from selenium.webdriver.firefox.options import Options

from covizu.utils.db_utils import *
from covizu.pangorider import *
from covizu.minimap2 import *

import subprocess
import argparse
import tempfile
import getpass
import json


def get_driver(download_folder, executable_path):
    """
    Instantiate remote control interface for Firefox web browser
    Debug code:print(driver.page_source)
    :param download_folder:  path to write downloaded files
    :param executable_path:  path to geckodriver executable
    :return:
    """
    profile = webdriver.FirefoxProfile()
    profile.set_preference('browser.download.folderList', 2)
    profile.set_preference('browser.download.manager.showWhenStarting', False)
    profile.set_preference("browser.download.dir", download_folder)
    profile.set_preference('browser.helperApps.alwaysAsk.force', False)
    profile.set_preference("browser.helperApps.neverAsk.saveToDisk", "application/octet-stream,octet-stream")
    profile.set_preference("browser.helperApps.neverAsk.openFile", "application/octet-stream,octet-stream")

    opts = Options()
    opts.headless = True  # opts.set_headless()
    assert opts.headless

    return webdriver.Firefox(firefox_profile=profile, options=opts, executable_path=executable_path)


def login(driver):
    """
    Use GISAID access credentials to login to database.
    :param driver: webdriver.Firefox object
    :return:
    """
    driver.get('https://www.epicov.org/epi3/cfrontend')
    time.sleep(15)  # seconds

    # read login credentials from Environment Variables
    try:
        user = os.environ['gisaid_u_variable']
        pw = os.environ['gisaid_pw_variable']
    except KeyError:
        # variables not set, get access credentials interactively
        user = getpass.getpass(prompt='GISAID username: ')
        pw = getpass.getpass(prompt='GISAID password: ')

    print('logging in')
    driver.execute_script('document.getElementById("elogin").value="{}"'.format(user))
    driver.execute_script('document.getElementById("epassword").value="{}"'.format(pw))
    time.sleep(5)

    # call javascript login function
    driver.execute_script('doLogin()')
    time.sleep(15)
    # navigate to corona virus page
    print('switching to Covid Homepage')
    element = driver.find_element_by_xpath("//*[contains(text(), 'EpiCoV™')]")
    element.click()
    time.sleep(15)

    print('navigating to CoV db')
    element = driver.find_element_by_xpath("//*[contains(text(), 'Browse')]")
    driver.execute_script("arguments[0].click();", element)
    #element.click() :TODO: this broke for some reason?
    time.sleep(5)

    return driver


def retrieve_genomes(driver, start, end, download_folder):
    """
    Retrieve genomes with a specified deposition date range in the GISAID database.
    Adding several time delays to avoid spamming the database.

    :param driver:  webdriver.Firefox object from login()
    :param start:  date in ISO format (yyyy-mm-dd)
    :param end:  date in ISO format (yyyy-mm-dd)
    :return:  path to file download
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

    # download seqs
    element = driver.find_element_by_xpath("//*[contains(text(), 'Download')]")
    driver.execute_script("arguments[0].click();", element)
    time.sleep(15)

    # switch to iframe to download
    driver.switch_to.frame(driver.find_element_by_tag_name("iframe"))
    print("Download")
    time.sleep(5)

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

    print('Downloading complete')

    downloaded_file = find_fasta(download_folder)
    while not downloaded_file.endswith( '.fasta'):
        time.sleep(15)
        downloaded_file = find_fasta(download_folder)
        print(downloaded_file)

    # reset browser
    time.sleep(30)
    driver.switch_to.default_content()
    element = driver.find_element_by_xpath("//button[@class='sys-event-hook sys-form-button']")
    element.click()

    #fix missing line breaks in-place
    retcode = subprocess.check_call(['sed', '-i', 's/([ACGT?])>hCo[Vv]/\1\\n>hCoV/g', downloaded_file])
    #fix spaces in header in-place
    retcode = subprocess.check_call(['sed', '-i', 's/ /_/g', downloaded_file])

    return downloaded_file


def find_fasta(download_folder):
    """ Get newest file in directory """
    downloaded_dir = os.listdir(download_folder)
    abs_dir = []
    for path in downloaded_dir:
        abs_dir.append(os.path.join(download_folder, path))
    downloaded_file = max(abs_dir, key=os.path.getmtime)
    return downloaded_file


def write_dbstats(db='data/gsaid.db'):
    """ write latest update string, with number of seqs """
    cur, conn = open_connection(db)
    numseqs = cur.execute('SELECT * FROM SEQUENCES')
    with open('data/dbstats.json', 'w') as jsonfile:
        data = {
            'lastupdate': date.today().isoformat(),
            'noseqs': len(numseqs.fetchall())
        }
        json.dump(data, jsonfile, indent=2)
    conn.close()


def parse_args():
    """ Command line interface """
    parser = argparse.ArgumentParser(
        description="Python3 Script to Automate retrieval of genomes deposited in a given day."
    )
    parser.add_argument(
        '--start', type=str, default=(date.today() - timedelta(days=1)).isoformat(),
        help="Start date to query database in ISO format (yyyy-mm-dd)."
    )
    parser.add_argument(
        '--end', type=str, default=date.today().isoformat(),
        help='End date to query database in ISO format (yyyy-mm-dd)'
    )
    parser.add_argument(
        '--destfile', type=str,
        default='data/gisaid-aligned.fa',
        help="Destination file to align and append downloaded sequences."
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
        '-r', '--ref', type= str, default='data/LR757995.fa',
        help='Path to reference fasta')
    parser.add_argument('--db', default = 'data/gsaid.db',
        help='Name of database.')
    parser.add_argument('--filterout', default = 'debug/filtered.log', type =argparse.FileType('w+'),
        help='Log for filtered seqs')
    parser.add_argument('--thread', default = 8,
        help='Number of threads for minimap')
    parser.add_argument('--minlen', help="<option> minimum sequence length, "
        "defaults to 29000nt.", default = 29000)
    site_packages = next(p for p in sys.path if 'site-packages' in p)
    parser.add_argument('--pangolindir', help ='PangoLEARN data dir, defaults to .../site_packages/panoLEARN/data/',
        default=site_packages+'/pangoLEARN/data/')
    parser.add_argument('--indicies', help='Indicies to keep, defaults to [265:29674]')
    return parser.parse_args()



if __name__ == '__main__':
    args = parse_args()
    indiciesToKeep = args.indicies

    driver = get_driver(download_folder=args.d.name, executable_path=args.binpath)
    driver = login(driver=driver)
    srcfile = retrieve_genomes(driver=driver, start=args.start, end=args.end,
                               download_folder=args.d.name)
    driver.quit()

    # align sequences in srcfile
    reflen = len(convert_fasta(open(args.ref))[0][1])
    print('Starting minimap')
    mm2 = minimap2(srcfile, ref=args.ref, nthread=args.thread,
                   minlen=args.minlen)

    # generate handle from mm2
    print('Reconstructing sequences')
    sequence_handle = stream_fasta(mm2, reflen=reflen)

    # insert aligned seqs into database
    iterate_handle(sequence_handle, args.db)

    # filter seqs here
    print('Filtering sequences for Lineage-typing')
    filtered_handle = filter_seqs(sequence_handle, args.filterout, max_prop_n=0.05, minlen=29000)
    print('Filtered {} seqs.'.format(str(len(filtered_handle))))

    # classify by PANGOLIN & insert processed seqs
    classify_and_insert(args.pangolindir+ 'decisionTreeHeaders_v1.joblib',
                        args.pangolindir+'decisionTree_v1.joblib',
                        filtered_handle, indiciesToKeep, args.db)

    print('Saving Lineages')
    insert_into_rawseqs(args.db, srcfile)
    write_dbstats()
