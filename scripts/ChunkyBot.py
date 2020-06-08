import os
import sys
import time
from datetime import date, timedelta, datetime

from selenium import webdriver
from selenium.webdriver.firefox.options import Options

import subprocess
import argparse
import tempfile
import getpass

from gotoh2 import *

todaystr = datetime.strftime(datetime.now(), '%Y-%m-%d')
cwd = os.getcwd()

def get_driver(download_folder, executable_path):
	"""
	Instantiate remote control interface for Firefox web browser
	
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
	time.sleep(5)

	#navigate to corona virus page
	print('navigating to CoV db')
	element = driver.find_element_by_xpath("//*[contains(text(), 'Browse')]")
	element.click()
	time.sleep(5)

	return driver

def find_prefix(driver):
	#find prefix variable
	element = driver.find_element_by_xpath("//div[@class='buttons container-slot']")
	return element.get_attribute('id').split('_')[1]

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

	# navigate to corona virus page
	print('navigating to CoV db')
	element = driver.find_element_by_xpath("//button[contains(text(), 'Reset')]")
	element.click()
	time.sleep(5)

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
	time.sleep(5)

	# switch to iframe to download
	driver.switch_to_frame(driver.find_element_by_tag_name("iframe"))
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
	downloaded_file = os.listdir(download_folder)[0]
	driver.switch_to.default_content()

	# reset browser
	element = driver.find_element_by_xpath("//button[contains(text(), 'Reset')]")
	element.click()
	return os.path.join(download_folder, downloaded_file)


def compare_dicts(old_fasta, new_fasta):
	old_dict = dict(convert_fasta(open(old_fasta)))
	new_dict = dict(convert_fasta(open(new_fasta)))
	block_diff={}
	for key, value in old_dict.items():
		try:
			newval = new_dict[key]
			if newval != value:
				block_diff[key] = newval
		except:
			block_diff[key] = value
	return block_diff

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
		help='Directory for Output file containing sequences that are modified or missing'
		)
	return parser.parse_args()


if __name__ == '__main__':
	args = parse_args()
	download_folder = args.d.name
	driver = get_driver(download_folder=download_folder, executable_path=args.binpath)
	driver = login(driver=driver)
	bfiles = ['data/baseline/'+ s for s in os.listdir(args.baselinedir)]
	difffasta ={}
	#predefined blocks one, two and three have <10k sequences
	srcfile = retrieve_genomes(driver=driver, start='2019-01-01', end='2020-04-16',
							   download_folder=download_folder)
	difffasta.update(compare_dicts(bfiles[0],srcfile))
	time.sleep(300)

	srcfile = retrieve_genomes(driver=driver, start='2020-04-17', end='2020-05-08',
							   download_folder=download_folder)
	difffasta.update(compare_dicts(bfiles[1],srcfile))
	time.sleep(300)

	srcfile = retrieve_genomes(driver=driver, start='2020-05-08', end='2020-05-15',
							   download_folder=download_folder)
	difffasta.update(compare_dicts(bfiles[2],srcfile))
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
		difffasta.update(compare_dicts(bfiles[osdircount],srcfile))
		osdircount +=1
		time.sleep(300)
		startdate = enddate +timedelta(1)
		#set up for next iteration

	#move latest file directly to baseline
	shutil.move(src, cwd+'/data/weeklydumps/baseline/'+blockname+ '.fasta')

	driver.quit()

	with open(args.outputdiffFasta, mode='w+') as out:
		for h,s in difffasta.items():
			out.write('>{}\n{}\n'.format(h, s))
