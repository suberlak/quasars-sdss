import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as  plt

fig,ax = plt.subplots(1,1) 

# k=3 , df = k-parameter = N-1 
df=3

x = np.linspace(chi2.ppf(0.0001,df),chi2.ppf(0.9999, df), 100)

ax.plot(np.log10(x), chi2.cdf(x, df),'r-', lw=5, alpha=0.6, label='chi2 cdf')
plt.xlim((-2,2))
plt.show()

