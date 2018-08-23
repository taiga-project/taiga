import sys
import numpy as np
import matplotlib.pyplot as pl
import scipy.stats as st

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

detector_RZ_angle = 38
R0 = 0.685
Z0 = 0.23
max_distance = 1e-4

results_folder = '../results/'+shotnumber+'_'+time+'/'+runnumber+'/'

try:
	R = np.loadtxt(results_folder+'rad.dat') #, delimiter="\t"
	Z = np.loadtxt(results_folder+'z.dat')
	T = np.loadtxt(results_folder+'tor.dat')
	
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

	plt.xlabel(r"$T \mathrm{[m]}$")
	plt.ylabel(r"$Z \mathrm{[m]}$")
	plt.savefig(results_folder+'detpy_'+shotnumber+'_'+time'.pdf')
	plt.clf()

except:
	print 'Invalid input folder: '+results_folder
