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
	R = np.loadtxt(results_folder+'rad.dat') #, delimiter="\t"
	Z = np.loadtxt(results_folder+'z.dat')
	T = np.loadtxt(results_folder+'tor.dat')
except:
	print 'Invalid input folder: '+results_folder

try:
	particle_on_detector = np.abs((Z-Z0) *np.tan((detector_RZ_angle/180.0)*np.pi) + (R-R0)) < max_distance
	
	x = T[particle_on_detector]
	y = Z[particle_on_detector]
	xmin, xmax = -4e-2, 4e-2
	ymin, ymax = 18e-2, 26e-2
	
	# Peform the kernel density estimate
	xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
	positions = np.vstack([xx.ravel(), yy.ravel()])
	values = np.vstack([x, y])
	kernel = st.gaussian_kde(values,0.05)
	f = np.reshape(kernel(positions).T, xx.shape)
	
	fig = pl.figure()
	ax = fig.gca()
	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)
	ax.set_aspect('equal')
	
	# Contourf plot
	#cfset = ax.contourf(xx, yy, f, 500, cmap='hot')  # Blues, afmhot; magma, inferno
	ax.hist2d(x, y, bins=250, range = [[xmin,xmax],[ymin,ymax]], cmap='afmhot')
	# Contour plot
	# Label plot

	pl.xlabel(r"$T \mathrm{ [m]}$")
	pl.ylabel(r"$Z \mathrm{ [m]}$")
except:
	print 'Unable to plot: '+results_folder

try:
	print 'Save plot to '+results_folder+'detpy_'+shotnumber+'_'+time+'.pdf'
	pl.savefig(results_folder+'detpy_'+shotnumber+'_'+time+'.pdf')
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
