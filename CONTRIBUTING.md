#The following Environment Variables need to be defined for GISAID downloads

gisaid_pw_variable='password'
gisaid_u_variable='username'

#The downloading scripts can be automated through crontab on linux:
0 0 * * * nohup python /home/covid/SeliniumAutobot.py >> /home/covid/Autobot.log 2>&1


