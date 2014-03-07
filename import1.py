# to take only RA, DEC, and errors on RA, Dec of stars from SDSS Stripe 82 , and
#  use them to find observations from CRTS for those stars
# http://nesssi.cacr.caltech.edu/DataRelease/multicone.html

# the problemn : CRTS takes only 100 objects per entry
# I have 1006849 rows of data, which would require over 10'000  manual entries! 
# there must be a better way of doing it 

# http://wwwsearch.sourceforge.net/mechanize/
# http://docs.python-requests.org/en/latest/user/quickstart/#passing-parameters-in-urls
# http://www.crummy.com/software/BeautifulSoup/
# http://stackoverflow.com/questions/5667699/python-urllib2-automatic-form-filling-and-retrieval-of-results
# https://www.distilled.net/blog/web-analytics/automating-search-query-downloads-from-webmaster-central/
# http://bcbio.wordpress.com/2010/01/02/automated-retrieval-of-expression-data-with-python-and-r/
# http://learnpythonthehardway.org/book/ex51.html

# try using 'requests'

import requests 
import numpy as np

arr=np.load('my_array1.npy')

ra=arr[:,0]   # all ra's 
dec=arr[:,1]  # all dec's

id1=np.arange(100)
rad1=np.zeros(100)+0.005  # search radius
ra1=arr[0:100,0] # first 100 ra's
dec1=arr[0:100,1] # first 100 dec's

np.savetxt('id_ra_dec_test.txt', np.column_stack((id1,ra1,dec1,rad1)),fmt='%11.4f')

