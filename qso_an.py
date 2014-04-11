# after performing all the averaging calculations on files located in the 
# directory with qso_var routine , we have 
#

#  make a list of output files 
# http://stackoverflow.com/questions/89228/calling-an-external-command-in-python
import os 
dire='QSO_try/'
cmd="cd "+dire
os.system(cmd)
os.system("pwd")
os.system("ls out* > out.list")
os.system("cd ..")
os.system("pwd")





