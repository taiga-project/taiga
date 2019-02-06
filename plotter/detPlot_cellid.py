import sys
import numpy as np
import matplotlib.pyplot as pl
import scipy.stats as st

max_distance = 1e-4

try:
	shotnumber = sys.argv[1];
except:
	print 'Shot number is not given as proper input parameter. Default value (17178) is used.'
	shotnumber = '17178'

try:
	time = sys.argv[2];
except:
	print 'Time is not given as proper input parameter. Default value (1005) is used.'
	time = '1005'

try:
	runnumber = sys.argv[3];
except:
	print 'Runnumber is not given as proper input parameter. Default value (402) is used.'
	runnumber = '402'

try:
	detector_par = sys.argv[4];
except:
	print 'Detector parameters is not given as proper input parameter. Default value (0.685,0.22,0,38,0) is used.'
	detector_par = '0.685,0.22,0,38,0'

detector_par_list = detector_par.split(',')

try:
	R0 = float(detector_par_list[0]);
except:
	R0 = 0.685

try:
	Z0 = float(detector_par_list[1]);
except:
	Z0 = 0.23

try:
	T0 = float(detector_par_list[2]);
except:
	T0 = 0.0

try:
	detector_RZ_angle = float(detector_par_list[3]);
except:
	detector_RZ_angle = 38.0

try:
	detector_ZT_angle = float(detector_par_list[4]);
except:
	detector_ZT_angle = 0.0
	
results_folder = 'results/'+shotnumber+'_'+time+'/'+runnumber+'/'

try:
	c = np.loadtxt(results_folder+'detcellcounter.dat')
except:
	print 'Invalid input folder: '+results_folder

C = c/c.max()

try:
	
	x = T[particle_on_detector]
	y = Z[particle_on_detector]
	xmin, xmax = -4e-2, 4e-2
	ymin, ymax = -4e-2, 4e-2
	
	fig = pl.figure()
	ax = fig.gca()
	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)
	ax.set_aspect('equal')
	
	# Contourf plot
	
	print C
	pl.imshow(C, cmap='afmhot')
	pl.xlabel(r"$T \mathrm{ [m]}$")
	pl.ylabel(r"$Z \mathrm{ [m]}$")
except:
	print 'Unable to plot: '+results_folder

try:
	print 'Save plot to '+results_folder+'detpy2_'+shotnumber+'_'+time+'.pdf'
	pl.savefig(results_folder+'detpy2_'+shotnumber+'_'+time+'.pdf')
	pl.clf()

except:
	print 'Unable to save to : '+results_folder
	pl.show()

xmin, xmax = -4e-2, 4e-2
ymin, ymax = -4e-2, 4e-2

#x = T[particle_on_detector]-T0
#y = (Z[particle_on_detector]-Z0) / np.sin( (detector_RZ_angle/180.0)*np.pi )
#nbins=250
#fig = pl.figure()
#ax = fig.gca()
#ax.set_xlim(xmin, xmax)
#ax.set_ylim(ymin, ymax)
#ax.set_aspect('equal')
#cfset = ax.hist2d(x, y, bins=nbins, range = [[xmin,xmax],[ymin,ymax]], cmap='afmhot')
#pl.colorbar(cfset[3]),ax=ax)
#pl.show()
