from selenium import webdriver
import time
from datetime import datetime, timedelta
import os
from selenium.webdriver.firefox.options import Options
import shutil 
import subprocess
from gotoh2 import *

todaystr = datetime.strftime(datetime.now(), '%Y-%m-%d')
cwd = os.getcwd()

print(todaystr + ' log')
print('===================')

destdir = cwd + '/data/weeklydumps/'+ todaystr +'/' 
bashcmd = 'mkdir -p ' + destdir
process = subprocess.Popen(bashcmd.split(), stdout=subprocess.PIPE)

#set profile settings (mostly so that it downloads automatically )
print('initializing options & profile')
profile = webdriver.FirefoxProfile()
profile.set_preference('browser.download.folderList', 2) 
profile.set_preference('browser.download.manager.showWhenStarting', False)
profile.set_preference("browser.download.dir", cwd + '/data/tmpdownloads')
profile.set_preference('browser.helperApps.alwaysAsk.force', False)
profile.set_preference("browser.helperApps.neverAsk.saveToDisk", "application/octet-stream,octet-stream")
profile.set_preference("browser.helperApps.neverAsk.openFile", "application/octet-stream,octet-stream")
opts = Options()
opts.set_headless()
assert opts.headless 

#initialize browser
print('launching browser ')
driver = webdriver.Firefox(firefox_profile=profile, options=opts, executable_path = '/usr/local/bin/geckodriver')
#print(driver)
driver.get('https://www.epicov.org/epi3/cfrontend')

time.sleep(5)
pw= os.environ['gisaid_pw_variable']
user= os.environ['gisaid_u_variable']
print('logging in')
driver.execute_script('document.getElementById("elogin").value="'+ user + '"')
driver.execute_script('document.getElementById("epassword").value="'+ pw + '"')

#call javascript login function 
driver.execute_script('doLogin()')
time.sleep(5)

#find prefix variable
element = driver.find_element_by_xpath("//div[@class='buttons container-slot']")
htmlid_as_list = element.get_attribute('id').split('_')
variable = htmlid_as_list[1]

#navigate to corona virus page
print('navigating to corona db')
element = driver.find_element_by_xpath("//*[contains(text(), 'Browse')]")
element.click()
time.sleep(5)

#change date range, and trigger selection change 
print('Selecting date')
yesterday = datetime.strftime(datetime.now() - timedelta(1), '%Y-%m-%d')

time.sleep(5)


def SelectAndDownload(start, end, filename):
	print('Selecting '+ start + ', '+ end)
	#takes start, end (two strings in format 2020-mm-dd) + filename and downloads sequences within that range 
	time_string =  '[id^="ce_' + variable + '"][id$="_input"]'
	#set start date 
	driver.execute_script("document.querySelectorAll('" + time_string +"')[2].value = '" + start +"'")
	driver.execute_script("document.querySelectorAll('" + time_string +"')[2].onchange()")
	#set end date 
	driver.execute_script("document.querySelectorAll('" + time_string +"')[3].value = '" + end +"'")
	driver.execute_script("document.querySelectorAll('" + time_string +"')[3].onchange()")
	driver.execute_script('''document.querySelectorAll('[id^="ce_''' + variable + '''"][id$=_input]')[2].onchange()''')
	#select all sequences 
	print('selecting all seqs')
	time.sleep(15)
	element = driver.find_element_by_xpath("//*[contains(text(), 'Total')]")
	count = element.get_attribute('innerHTML').split()[1].replace(',','')
	if int(count) > 10000:
		time.sleep(15)  
	checkbox = driver.find_element_by_xpath("//span[@class='yui-dt-label']/input[@type='checkbox']")
	checkbox.click()
	time.sleep(5)
	#click download button, look for button that contains Download txt
	element = driver.find_element_by_xpath("//*[contains(text(), 'Download')]")
	driver.execute_script("arguments[0].click();", element)
	#element.click()
	time.sleep(5)
	#switch to iframe to download 
	driver.switch_to_frame(driver.find_element_by_tag_name("iframe"))
	# download
	print("Download")
	time.sleep(5)
	button = driver.find_element_by_xpath("//*[contains(text(), 'Download')]//ancestor::div[@style='float: right']")
	button.click()
	time.sleep(5)
	# wait for download to complete
	downloading = True
	while downloading:
		if os.listdir(cwd +'/data/tmpdownloads/') == []:
			time.sleep(10)
		for file in os.listdir(cwd + '/data/tmpdownloads/'):
			if file.endswith('.part'):
				time.sleep(5)
				downloading = True
		if len(os.listdir(cwd+'/data/tmpdownloads/')) == 1:
			downloading = False

	print('Downloading complete, moving files')
	#move file
	yesterday = datetime.strftime(datetime.now() - timedelta(1), '%Y-%m-%d')
	downloadedfile = os.listdir(cwd+'/data/tmpdownloads/')[0]
	source = cwd+'/data/tmpdownloads/' + downloadedfile
	destination = destdir + 'GISAID-' +  filename + '.fasta'
	shutil.move(source, destination)
	driver.switch_to.default_content()
	element = driver.find_element_by_xpath("//button[@class='sys-event-hook sys-form-button']")
	element.click()

#predefined blocks one, two and three have <10k sequences 
SelectAndDownload('2019-01-01', '2020-04-16', '0101_0416')
time.sleep(300)
driver.switch_to.default_content()
SelectAndDownload('2020-04-17', '2020-05-08', '0417_0508')
time.sleep(300)
driver.switch_to.default_content()
SelectAndDownload('2020-05-08', '2020-05-15', '0508_0515')
time.sleep(300)
driver.switch_to.default_content()
#next do the chunks afterwards in weekly increments 
#find how many new chunks we need to get 
startdate = datetime(2020,5,16)
today = datetime.now()
delta = (today- startdate).days
#find number of weekly blocks
newchunks = delta//7 + 1 

for count in range(1, newchunks + 1):
	#get the chunk dates
	enddate = startdate + timedelta(6)
	startmonthstr = str(startdate.month) if startdate.month > 9 else '0'+ str(startdate.month)
	startdaystr = str(startdate.day) if startdate.day > 9 else '0' + str(startdate.day)
	endmonthstr = str(enddate.month) if enddate.month > 9 else '0'+ str(enddate.month)
	enddaystr = str(enddate.day) if enddate.day > 9 else '0' + str(enddate.day)
	blockname = startmonthstr+startdaystr + '_' + endmonthstr+enddaystr
	SelectAndDownload(datetime.strftime(startdate, '%Y-%m-%d'), datetime.strftime(enddate, '%Y-%m-%d'), blockname)
	time.sleep(20)
	driver.switch_to.default_content()
	startdate = enddate +timedelta(1)
	#set up for next iteration 


driver.quit()

#move latest file directly to baseline

shutil.move(destdir + 'GISAID' +blockname+ '.fasta', cwd+'/data/weeklydumps/baseline')

#Load and compare Fasta seqs

all_diff_headers =[]

#compare blocks:

def import_dicts(old_fasta, new_fasta):
	base_fasta = convert_fasta(open(old_fasta))
	new_fasta = convert_fasta(open(new_fasta))
	return (dict(base_fasta), dict(new_fasta))

def compare_dicts(old_dict, new_dict):
	block_diff=[]
	for key, value in old_dict.iteritems():
		try:
			newval = new_dict[key]
			if newval != value:
				block_diff.append(key)
		except:
			block_diff.append(key)
	return block_diff

for file in os.listdir(cwd +'/data/weeklydumps/baseline'):
	old_dict, new_dict = import_dicts(cwd +'/data/weeklydumps/baseline/'+file, destdir +file)
	for diff_header in compare_dicts(old_dict, new_dict):
		all_diff_headers.append(diff_header)

#write report
reportFile = open(cwd+'/debug/'+ todaystr + '_report.txt', 'w+')

for header in all_diff_headers:
	reportFile.write(header +'\n')

#destroy additional files in no differences reported
if len(all_diff_headers) == 0:
	reportFile.write('No Differences Found')
	shutil.rmtree(destdir)
 
reportFile.close()



