import sys
import os
import numpy
from scipy.interpolate import interp1d
import h5py
import matplotlib.pyplot as plt

numberOfCmdArgs = len(sys.argv)

rootFolder = os.path.abspath(os.getcwd())
print(rootFolder)
resultFolder = rootFolder+'/results';
fieldFolder  = rootFolder+'/input/fieldGrid';
cdbFolder    = rootFolder+'/input/cdb';
ionProfileFolder =rootFolder+'/input/ionProf';

print(numberOfCmdArgs)
if (numberOfCmdArgs > 1):
    shotAndTime = sys.argv[1]
else:
    shotAndTime = '11774_1000'

if (numberOfCmdArgs > 2):
    runnumber = sys.argv[2]
else:
    runnumber = '0';

workingFolder = resultFolder+'/'+shotAndTime+'/'+runnumber

#ionisation
ionR = numpy.genfromtxt(fname=ionProfileFolder+'/'+shotAndTime+'/rad.dat')
ionY = numpy.genfromtxt(fname=ionProfileFolder+'/'+shotAndTime+'/ionyeald.dat')
ionX1 = numpy.gradient(ionY);
ionF = interp1d(ionR, -ionX1/numpy.amax(numpy.absolute(ionX1)));

plt.plot(ionR, ionF(ionR))

try:
    print ('Save plot to '+workingFolder+'/ion_'+shotAndTime+'.pdf')
    plt.savefig(workingFolder+'/ion_'+shotAndTime+'.pdf')
    plt.clf()

except:
    print ('Unable to save to : '+workingFolder+'/ion_'+shotAndTime+'.pdf')
    plt.show()

# trajectories
print(resultFolder+'/'+shotAndTime+'/'+runnumber+'/t_rad.dat')
t_rad = numpy.genfromtxt(fname=workingFolder+'/t_rad.dat')
t_z   = numpy.genfromtxt(fname=workingFolder+'/t_z.dat')
t_tor = numpy.genfromtxt(fname=workingFolder+'/t_tor.dat')

numberOfIons = numpy.size(t_rad, 1)

print(fieldFolder+'/'+shotAndTime+'/rcord.dat')
rcord = numpy.genfromtxt(fname=fieldFolder+'/'+shotAndTime+'/rcord.dat')
zcord = numpy.genfromtxt(fname=fieldFolder+'/'+shotAndTime+'/zcord.dat')
psi2  = numpy.genfromtxt(fname=fieldFolder+'/'+shotAndTime+'/psi2.dat')

shotnumber = shotAndTime.split('_')[0]
timeSlice  = shotAndTime.split('_')[1]

hdf5Path = cdbFolder+'/'+shotnumber+'/EFITXX/EFITXX.1.h5'
print(hdf5Path)
hdf5File = h5py.File(hdf5Path, 'r')
hdf5Times = hdf5File.get('time').value
hdf5PsiBoundary = hdf5File.get('output/globalParameters/psiBoundary').value
hdf5PsiAxis = hdf5File.get('output/globalParameters/psiAxis').value
hdf5BoundaryR = hdf5File.get('output/separatrixGeometry/boundaryCoordsR').value
hdf5BoundaryZ = hdf5File.get('output/separatrixGeometry/boundaryCoordsZ').value

hdf5Index = int(numpy.argwhere(hdf5Times == float(timeSlice)/1000.))

psiAxis = hdf5PsiAxis[hdf5Index]
psiBoundary = hdf5PsiBoundary[hdf5Index]
boundaryR = hdf5BoundaryR[hdf5Index,:]
boundaryZ = hdf5BoundaryZ[hdf5Index,:]

PSI = psi2/psiAxis

R, Z =  numpy.meshgrid(rcord, zcord)

if (numberOfCmdArgs > 3):
    t_index = sys.argv[2]
    t_rad = t_rad[:,t_index]
    t_z   = t_z  [:,t_index]
    t_tor = t_tor[:,t_index]

#ion histogram (real)
fig = plt.figure()
plt.hist(t_rad[0,:], density=True, bins=50)
plt.xlabel(r'$R$ [m]')
plt.title (r'COMPASS #'+shotnumber+'\n ($t = $'+timeSlice+' ms)')

try:
    print ('Save plot to '+workingFolder+'/start_'+shotAndTime+'.pdf')
    plt.savefig(workingFolder+'/start_'+shotAndTime+'.pdf')
    plt.clf()

except:
    print ('Unable to save to : '+workingFolder+'/start_'+shotAndTime+'.pdf')
    plt.show()

#detector
plt.clf()
fig, axs = plt.subplots(2, 2)
axs[0,0].plot(t_rad[-1,:], t_z[-1,:], '.')
#axs[0,0].set_xlabel(r'$R$ [m]')
axs[0,0].set_ylabel(r'$Z$ [m]')
axs[0,0].yaxis.set_ticks_position('both')

axs[0,1].plot(t_tor[-1,:], t_z[-1,:], '.')
axs[0,1].set_xlabel(r'$T$ [m]')
axs[0,1].set_ylabel(r'$Z$ [m]')
axs[0,1].yaxis.tick_right()
axs[0,1].yaxis.set_ticks_position('both')

axs[1,0].plot(t_rad[-1,:], t_tor[-1,:], '.')
axs[1,0].set_xlabel(r'$R$ [m]')
axs[1,0].set_ylabel(r'$T$ [m]')
axs[1,0].yaxis.set_ticks_position('both')

axs[1,1].axis('off')

fig.suptitle(r'COMPASS #'+shotnumber+' ($t = $'+timeSlice+' ms)')
#plt.colorbar()

try:
    print ('Save plot to '+workingFolder+'/end_'+shotAndTime+'.pdf')
    plt.savefig(workingFolder+'/end_'+shotAndTime+'.pdf')
    plt.clf()

except:
    print ('Unable to save to : '+workingFolder+'/end_'+shotAndTime+'.pdf')
    plt.show()

#2D figure
plt.clf()
fig = plt.figure()
#plt.contourf(R, Z, -PSI, levels=numpy.arange(-1.05,2,.05), cmap='Spectral')
#plt.contourf(R, Z, PSI, levels=numpy.arange(-2,2,.05), cmap='coolwarm')
#plt.contourf(R, Z, PSI, levels=numpy.arange(-2,2,.05), cmap='twilight_shifted')
plt.contourf(R, Z, PSI, levels=numpy.arange(-2,2,.1), cmap='GnBu')
plt.plot(boundaryR, boundaryZ, color=(0,0,0))
for i in range(numberOfIons):
    i1 = 1
    i2 = ionF(t_rad[1,i])
    plt.plot([rcord[-1], t_rad[1,i]], [t_z[0,i], t_z[0,i]], color=(i1, i1, 0))
    if i2<0:
        i2=0
    if i2>1:
        i2=1
    plt.plot(t_rad[:,i], t_z[:,i], color=(i2, 0, 0))

plt.xlabel(r'$R$ [m]')
plt.ylabel(r'$Z$ [m]')
plt.title (r'COMPASS #'+shotnumber+'\n ($t = $'+timeSlice+' ms)')
ax = fig.add_subplot(111)
ax.set_aspect('equal', 'box')
#plt.colorbar()

try:
    print ('Save plot to '+workingFolder+'/traj_'+shotAndTime+'.pdf')
    plt.savefig(workingFolder+'/traj_'+shotAndTime+'.pdf')
    plt.clf()

except:
    print ('Unable to save to : '+workingFolder+'/traj_'+shotAndTime+'.pdf')
    plt.show()

