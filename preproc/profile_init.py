import numpy as np
import h5py
import sys
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


def prof_psi(x,a,b,c):
	return  a * np.exp(-(np.sqrt(x)/b)**c)

def prof_r(r,a,b,c):
	return  a * np.exp(-(r/b)**c)

def prof_r2(r,a,b,c):
	return  a * np.exp(-(r/b)**3.0)+c

def prof_r3(r,b,c):
	return np.exp(-(r/b)**c)

shotnumber =  sys.argv[1];

f=h5py.File('../input/cdb/'+shotnumber+'/THOMSON/ne.1.h5','r')
ne = f['ne']
f=h5py.File('../input/cdb/'+shotnumber+'/THOMSON/Te.1.h5','r')
Te = f['Te']

f=h5py.File('../input/cdb/'+shotnumber+'/THOMSON/ne_err.1.h5','r')
ne_err = f['ne_err']
f=h5py.File('../input/cdb/'+shotnumber+'/THOMSON/Te_err.1.h5','r')
Te_err = f['Te_err']

f=h5py.File('../input/cdb/'+shotnumber+'/THOMSON/psi_n.1.h5','r')
psi = f['psi_n']
f=h5py.File('../input/cdb/'+shotnumber+'/THOMSON/TS_record_time.1.h5','r')
t = f['TS_record_time']
t_int = (np.rint(t)).astype(int)

#np.savetxt('../input/tsProf/'+shotnumber+'.time',t_int,fmt='%d')

f=h5py.File('../input/cdb/'+shotnumber+'/EFITXX/EFITXX.1.h5','r')
psi_efit = f['output/radialProfiles/normalizedPoloidalFlux']
r_efit = f['output/radialProfiles/r']

t_all=[]

for i in range(1, len(t)):
	try:
		efit_min_index = np.argmin(psi_efit[i])
		psi_to_r = interp1d(psi_efit[i][efit_min_index:], r_efit[i][efit_min_index:],  kind='cubic')

		poly_psi_to_r = np.polyfit(psi_efit[i], r_efit[i], 4)
		poly_r_to_psi = np.polyfit(r_efit[i], psi_efit[i], 3)

		ts_min_index = np.nanargmin(psi[i])#np.nanargmax(ne[i])

		max_psi_limit = 0.9 #0.9
		valid_ts = np.logical_and(~np.isnan(ne[i]) ,  psi[i] < max_psi_limit);
		#valid_ts = psi[i][valid_ts] < max_psi_limit

		ne_ts = ne[i][valid_ts][ts_min_index:];
		Te_ts = Te[i][valid_ts][ts_min_index:];
		psi_ts = psi[i][valid_ts][ts_min_index:];	
		r_ts = psi_to_r(psi_ts)

		ne_filter = psi_ts>0.1

		#ne_norm = ne_ts/ne_ts[0]
		#Te_norm = Te_ts/Te_ts[0]
		#r_ne_half = np.interp(0.5,ne_norm[ne_filter],r_ts[ne_filter])
		#r_Te_half = np.interp(0.5,Te_norm[ne_filter],r_ts[ne_filter])

		ne_coeff, ne_cc = curve_fit(prof_r3,r_ts[ne_filter],ne_ts[ne_filter]/ne_ts[ne_filter][0], p0=(0.7,4.))
		Te_coeff, Te_cc = curve_fit(prof_r3,r_ts[ne_filter],Te_ts[ne_filter]/Te_ts[ne_filter][0], p0=(0.7,4.))
		r_minmax = np.polyval(poly_psi_to_r,np.array([0.0, 2.0]))
		r_more=np.linspace(r_minmax[0],0.8,200)
		ne_more = prof_r3(r_more,ne_coeff[0],ne_coeff[1])*ne_ts[ne_filter][0]
		Te_more = prof_r3(r_more,Te_coeff[0],Te_coeff[1])*Te_ts[ne_filter][0]

		psi_more = np.polyval(poly_r_to_psi,r_more)

		if False:	
			ne_more0 = prof_r3(r_more,ne_coeff[0],ne_coeff[1])
			Te_more0 = prof_r3(r_more,Te_coeff[0],Te_coeff[1])
			plt.figure(1)
			plt.plot(r_ts,ne_ts/ne_ts[ne_filter][0],'ro')
			plt.plot(r_ts[ne_filter],ne_ts[ne_filter]/ne_ts[ne_filter][0],'x')
			plt.plot(r_more,ne_more0,'r-')
			plt.plot(r_ts,Te_ts/Te_ts[ne_filter][0],'bo')
			plt.plot(r_ts[ne_filter],Te_ts[ne_filter]/Te_ts[ne_filter][0],'x')
			plt.plot(r_more,Te_more0,'b-')
			plt.figure(2)
			plt.plot(r_ts,ne_ts,'ro')
			plt.plot(r_more,ne_more,'r-')
			plt.figure(3)
			plt.plot(r_ts,Te_ts,'ro')
			plt.plot(r_more,Te_more,'r-')
			plt.show()

		if False:	
			plt.plot(psi_more,ne_more,'r-')
			plt.errorbar(psi[i],ne[i], yerr=ne_err[i]*3, fmt='ro')
			plt.show()

		np.savetxt('../input/tsProf/ne_'+shotnumber+'_'+str(t_int[i])+'.dat',ne[i])
		np.savetxt('../input/tsProf/Te_'+shotnumber+'_'+str(t_int[i])+'.dat',Te[i])
		o = np.ones_like(psi_more)
		renate_iio_data = np.matrix([psi_more,Te_more/1e3,ne_more/1e19,2.0*o,5.0*o])
		np.savetxt('../input/tsProf/nt_compass'+shotnumber+'_'+str(t_int[i])+'.txt',renate_iio_data.getT())
		t_all.append(t_int[i])
	except:
		print 'Error at '+str(t_int[i]), sys.exc_info()[0]

np.savetxt('../input/tsProf/'+shotnumber+'.time',t_all,fmt='%d')

