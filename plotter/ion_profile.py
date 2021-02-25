import sys
import numpy
from scipy.interpolate import interp1d
import h5py
import matplotlib.pyplot as plt

numberOfCmdArgs = len(sys.argv)

rootFolder = '..'
rootFolder = '/home/matyi/work/taiga_local'
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

ionR = numpy.genfromtxt(fname=ionProfileFolder+'/'+shotAndTime+'/rad.dat')
ionY = numpy.genfromtxt(fname=ionProfileFolder+'/'+shotAndTime+'/ionyeald.dat')
ionX1 = numpy.gradient(ionY);
ionF = interp1d(ionR, -ionX1/numpy.amax(numpy.absolute(ionX1)));

numberOfIons = numpy.size(t_rad, 1)
plt.plot(ionR, ionF(ionR))

try:
    print ('Save plot to '+workingFolder+'/ion_'+shotAndTime+'.pdf')
    plt.savefig(workingFolder+'/ion_'+shotAndTime+'.pdf')
    plt.clf()

except:
    print ('Unable to save to : '+workingFolder+'/ion_'+shotAndTime+'.pdf')
    plt.show()

