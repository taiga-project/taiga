import numpy as np
import h5py
import sys
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


def l_mode_profile(r,b,c):
	return np.exp(-(r/b)**c)

try:
	shotnumber =  sys.argv[1];
except:
	print ('Shot number is not given as proper input parameter. Default value (11774) is used.'
	shotnumber = '11774'

try:
	f=h5py.File('input/cdb/'+shotnumber+'/THOMSON/ne.1.h5','r')
	ne = f['ne']
	f=h5py.File('input/cdb/'+shotnumber+'/THOMSON/Te.1.h5','r')
	Te = f['Te']

	f=h5py.File('input/cdb/'+shotnumber+'/THOMSON/ne_err.1.h5','r')
	ne_err = f['ne_err']
	f=h5py.File('input/cdb/'+shotnumber+'/THOMSON/Te_err.1.h5','r')
	Te_err = f['Te_err']

	f=h5py.File('input/cdb/'+shotnumber+'/THOMSON/psi_n.1.h5','r')
	psi = f['psi_n']
	f=h5py.File('input/cdb/'+shotnumber+'/THOMSON/TS_record_time.1.h5','r')
	t = f['TS_record_time']
	t_int = (np.rint(t)).astype(int)

	#np.savetxt('input/tsProf/'+shotnumber+'.time',t_int,fmt='%d')

	f=h5py.File('input/cdb/'+shotnumber+'/EFITXX/EFITXX.1.h5','r')
	psi_efit = f['output/radialProfiles/normalizedPoloidalFlux']
	r_efit = f['output/radialProfiles/r']

	t_all=[]

	for i in range(1, len(t)):
		try:
			efit_min_index = np.argmin(psi_efit[i])
			psi_to_r = interp1d(psi_efit[i][efit_min_index:], r_efit[i][efit_min_index:],  kind='cubic')

			poly_psi_to_r = np.polyfit(psi_efit[i], r_efit[i], 4)
			poly_r_to_psi = np.polyfit(r_efit[i], psi_efit[i], 3)

			ts_min_index = np.nanargmin(psi[i])

			max_psi_limit = 0.95 #0.9
			valid_ts = np.logical_and(~np.isnan(ne[i]) ,  psi[i] < max_psi_limit);

			ne_ts = ne[i][valid_ts][ts_min_index:];
			Te_ts = Te[i][valid_ts][ts_min_index:];
			psi_ts = psi[i][valid_ts][ts_min_index:];	
			r_ts = psi_to_r(psi_ts)

			ne_filter = psi_ts>0.05 #0.1

			ne_coeff, ne_cc = curve_fit(l_mode_profile,r_ts[ne_filter],ne_ts[ne_filter]/ne_ts[ne_filter][0], p0=(0.7,4.))
			Te_coeff, Te_cc = curve_fit(l_mode_profile,r_ts[ne_filter],Te_ts[ne_filter]/Te_ts[ne_filter][0], p0=(0.7,4.))
			r_minmax = np.polyval(poly_psi_to_r,np.array([0.0, 2.0]))
			r_more=np.linspace(r_minmax[0],0.8,200)
			ne_more = l_mode_profile(r_more,ne_coeff[0],ne_coeff[1])*ne_ts[ne_filter][0]
			Te_more = l_mode_profile(r_more,Te_coeff[0],Te_coeff[1])*Te_ts[ne_filter][0]

			psi_more = np.polyval(poly_r_to_psi,r_more)

			if False:	
				ne_more0 = l_mode_profile(r_more,ne_coeff[0],ne_coeff[1])
				Te_more0 = l_mode_profile(r_more,Te_coeff[0],Te_coeff[1])
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
				#plt.show()
				plt.clf()

			if True:
				plt.plot(r_ts,ne_ts,'bo')
				plt.plot(r_more,ne_more,'r-')
				plt.xlim((0,2))
				plt.ylim((0,1e20))
				plt.xlabel(r"$\rho_{norm}$")
				plt.ylabel(r"$n_e \mathrm{[m^{-3}]}$")
				plt.savefig('input/renate110/ne_r_'+shotnumber+'_'+str(t_int[i])+'.pdf')
				plt.clf()
				
				plt.plot(r_ts,Te_ts,'bo')
				plt.plot(r_more,Te_more,'r-')
				plt.xlim((0,2))
				plt.ylim((0,2e3))
				plt.xlabel(r"$\rho_{norm}$")
				plt.ylabel(r"$T_e \mathrm{[eV]}$")
				plt.savefig('input/renate110/te_r_'+shotnumber+'_'+str(t_int[i])+'.pdf')
				plt.clf()

				plt.plot(psi_more,ne_more,'r-')
				plt.fill_between(psi[i], ne[i]-ne_err[i]*3, ne[i]+ne_err[i]*3,facecolor='lightgrey')
				plt.errorbar(psi[i],ne[i], yerr=ne_err[i], fmt='bo-')
				plt.xlim((0,2))
				plt.ylim((0,1e20))
				plt.xlabel(r"$\psi$")
				plt.ylabel(r"$n_e \mathrm{[m^{-3}]}$")
				plt.savefig('input/renate110/ne_'+shotnumber+'_'+str(t_int[i])+'.pdf')
				plt.clf()

				plt.plot(psi_more,Te_more,'r-')
				plt.fill_between(psi[i], Te[i]-Te_err[i]*3, Te[i]+Te_err[i]*3,facecolor='lightgrey')
				plt.errorbar(psi[i],Te[i], yerr=Te_err[i], fmt='bo-')
				plt.xlim((0,2))
				plt.ylim((0,2e3))
				plt.xlabel(r"$\psi$")
				plt.ylabel(r"$T_e \mathrm{[eV]}$")
				plt.savefig('input/renate110/te_'+shotnumber+'_'+str(t_int[i])+'.pdf')
				plt.clf()

			np.savetxt('input/tsProf/ne_'+shotnumber+'_'+str(t_int[i])+'.dat',ne[i])
			np.savetxt('input/tsProf/Te_'+shotnumber+'_'+str(t_int[i])+'.dat',Te[i])
			o = np.ones_like(psi_more)
			renate_iio_data = np.matrix([psi_more,Te_more/1e3,ne_more/1e19,2.0*o,5.0*o])
			np.savetxt('input/renate110/nt_compass'+shotnumber+'_'+str(t_int[i])+'.txt',renate_iio_data.getT())
			t_all.append(t_int[i])
		except:
			print 'There is no valid Thomson profile at t = '+str(t_int[i])+' s'

	np.savetxt('input/tsProf/'+shotnumber+'.time',t_all,fmt='%d')
except:
	print 'There is no Thomson data in '+'input/cdb/'+shotnumber+'/THOMSON/'
