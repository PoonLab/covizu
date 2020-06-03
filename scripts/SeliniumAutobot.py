from selenium import webdriver
import time
from datetime import datetime, timedelta
import os
from selenium.webdriver.firefox.options import Options
import shutil  
import subprocess


cwd = os.getcwd()
download_folder = cwd+ '/data/tmpdownloads/'
today = datetime.strftime(datetime.now(), '%Y-%m-%d')

print(today + ' log')
print('===================')

#set profile settings (mostly so that it downloads automatically )
print('initializing options & profile')
profile = webdriver.FirefoxProfile()
profile.set_preference('browser.download.folderList', 2) 
profile.set_preference('browser.download.manager.showWhenStarting', False)
profile.set_preference("browser.download.dir", download_folder)
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

time.sleep(15)
#Reads login credentials from Environment Variables
pw= os.environ['gisaid_pw_variable']
user= os.environ['gisaid_u_variable']
print('logging in')
driver.execute_script('document.getElementById("elogin").value="'+user+ '"')
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


time_string =  '[id^="ce_' + variable + '"][id$="_input"]'
driver.execute_script("document.querySelectorAll('" + time_string +"')[2].value = '" + yesterday +"'")
driver.execute_script("document.querySelectorAll('" + time_string +"')[2].onchange()")
driver.execute_script('''document.querySelectorAll('[id^="ce_''' + variable + '''"][id$=_input]')[2].onchange()''')

time.sleep(5)


#select all sequences 
print('selecting all seqs')
checkbox = driver.find_element_by_xpath("//span[@class='yui-dt-label']/input[@type='checkbox']")
checkbox.click()
time.sleep(5)

#download seqs
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
	for file in os.listdir(download_folder):
		downloading = False
		if file.endswith('.part'):
			time.sleep(5)
			downloading = True

#move file
print('Downloading complete, moving files')

yesterday = datetime.strftime(datetime.now() - timedelta(1), '%Y-%m-%d')

downloadedfile = os.listdir(download_folder)[0]
source = download_folder + downloadedfile
destination = cwd + '/data/GISAID-' +  datetime.strftime(datetime.now(), '%Y-%m-%d') + '.fasta'

#call Sed command to fix  missing line characters

outfile = open(destination, 'w')
bashcmd = "sed s/[ACGT?]>hCoV/\\n>hCov/g " + source
process = subprocess.Popen(bashcmd.split(), stdout=outfile)
outfile.close()
driver.quit()

#run update script 
print('Running updater script')
bashcmd = '/usr/bin/python ' + cwd + '/scripts/update.py -ref '+ cwd + '/data/NC_045512.fa '+ destination+' ' + cwd +'/data/gisaid-aligned.fa'
#print(bashcmd)
process = subprocess.Popen(bashcmd.split(), stdout=subprocess.PIPE)
output, error = process.communicate()
print('Updater output')
print('===========')
print(output +'\n')

#write latest update string
file= open(cwd + '/data/lastupdate.json', 'w')
file.write('var lastupdate="' + datetime.strftime(datetime.now(), '%Y-%m-%d') + '\";')
file.close()
