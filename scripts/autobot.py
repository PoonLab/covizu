import os
import sys
import time
from datetime import date, timedelta

from selenium import webdriver
from selenium.webdriver.firefox.options import Options

import subprocess
import argparse
import tempfile
import getpass


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

        # navigate to corona virus page
	print('navigating to CoV db')
	element = driver.find_element_by_xpath("//*[contains(text(), 'Browse')]")
	element.click()
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

	# reset browser
	time.sleep(30)
	driver.switch_to.default_content()
	element = driver.find_element_by_xpath("//button[@class='sys-event-hook sys-form-button']")
	element.click()

	return os.path.join(download_folder, downloaded_file)


def update_local(srcfile, destfile):
	"""
	Call update.py for pairwise alignment and appending of new sequences to local file.
	:param srcfile:  path to FASTA file with downloaded genomes
	:param destfile:  path to local FASTA file of aligned genomes
	:return:
	"""
	# fix missing line breaks in-place
	retcode = subprocess.check_call(['sed', '-i', 's/([ACGT?])>hCo[Vv]/\1\\n>hCoV/g', srcfile])

	# call updater script
	process = subprocess.Popen(
		[sys.executable, 'scripts/update.py', srcfile, destfile],
		stdout=subprocess.PIPE
	)
	output, error = process.communicate()
	print('Updater output')
	print('==============')
	print(output + b'\n')

	# write latest update string
	with open('data/lastupdate.json', 'w') as jsonfile:
		jsonfile.write('var lastupdate="{}";'.format(date.today().isoformat()))


def parse_args():
	""" Command line interface """
	parser = argparse.ArgumentParser(
		description="Python3 Script to Automate retrieval of genomes deposited in a given day."
	)
	parser.add_argument(
		'start', type=str, default=(date.today() - timedelta(days=1)).isoformat(),
		help="Start date to query database in ISO format (yyyy-mm-dd)."
	)
	parser.add_argument(
		'end', type=str, default=date.today().isoformat(),
		help='End date to query database in ISO format (yyyy-mm-dd)'
	)
	parser.add_argument(
		'destfile', type=str,
		default='data/gisaid-aligned.fa',
		help="Destination file to align and append downloaded sequences."
	)
	parser.add_argument(
		'-d', '-dir', type=str, default=tempfile.TemporaryDirectory().name,
		help="Temporary directory to download files."
	)
	parser.add_argument(
		'-b', '--binpath', type=str, default='/usr/local/bin/geckodriver',
		help='Path to geckodriver binary executable'
	)
	return parser.parse_args()


if __name__ == '__main__':
	args = parse_args()
	driver = get_driver(download_folder=args.d, executable_path=args.binpath)
	driver = login(driver=driver)
	srcfile = retrieve_genomes(driver=driver, start=args.start, end=args.end,
							   download_folder=args.d)
	update_local(srcfile=srcfile, destfile=args.destfile)
	driver.quit()
