import numpy
import matplotlib.pyplot as plt


levels=numpy.linspace(-0.2, 0.2, 20)
levels2=numpy.linspace(-2, 2, 20)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

x = numpy.genfromtxt("exported_fieldR.dat", dtype=None)
R = x[:,0]
Z = x[:,1]
BR = x[:,2]
ax1.tricontourf(R, Z, BR, levels=levels)

x = numpy.genfromtxt("exported_fieldZ.dat", dtype=None)
R = x[:,0]
Z = x[:,1]
BZ = x[:,2]
ax2.tricontourf(R, Z, BZ, levels=levels)

x = numpy.genfromtxt("exported_fieldT.dat", dtype=None)
R = x[:,0]
Z = x[:,1]
BT = x[:,2]
ax3.tricontourf(R, Z, BT, levels=levels2)

Brad = numpy.genfromtxt("input/fieldGrid/17178_1097/brad.dat", dtype=None)
Btor = numpy.genfromtxt("input/fieldGrid/17178_1097/btor.dat", dtype=None)
R = numpy.genfromtxt("input/fieldGrid/17178_1097/rcord.dat", dtype=None)
Z = numpy.genfromtxt("input/fieldGrid/17178_1097/zcord.dat", dtype=None)
ax4.contourf(R, Z, Btor.T, levels=levels2)

plt.show()

