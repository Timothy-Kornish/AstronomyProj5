import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pickle

mags_b = pickle.load(open('/Volumes/TIM/Lab5/trans/mag_b.txt','rb'))
time_b = pickle.load(open('/Volumes/TIM/Lab5/trans/time_b.txt', 'rb'))
mags_a = pickle.load(open('/Volumes/TIM/Lab5/trans/mag_a.txt','rb'))
time_a = pickle.load(open('/Volumes/TIM/Lab5/trans/time_a.txt', 'rb'))

plt.scatter(time_b,mags_b)
plt.scatter( time_a,mags_a)
plt.xlabel('UT (hours)')   
plt.ylabel('Flux (Counts/Second)')
plt.xlim(1.4,7.4)
plt.ylim(18000,40000)
plt.show()