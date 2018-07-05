import numpy as np
import h5py
import sys
from scipy.optimize import curve_fit

xdata = np.array([-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9])
ydata = np.array([0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001])

def func(x, p1,p2):
  return p1*np.cos(p2*x) + p2*np.sin(p1*x)

popt, pcov = curve_fit(func, xdata, ydata,p0=(1.0,0.2))
print popt
#raw_input("Press any Key to Exit")

shotnumber =  sys.argv[1];

f=h5py.File('../input/cdb/'+shotnumber+'/THOMSON/ne.1.h5','r')
ne = f['ne']
f=h5py.File('../input/cdb/'+shotnumber+'/THOMSON/Te.1.h5','r')
Te = f['Te']
f=h5py.File('../input/cdb/'+shotnumber+'/THOMSON/TS_record_time.1.h5','r')
t = f['TS_record_time']
