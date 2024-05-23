import sys
import os
import numpy
from scipy.interpolate import interp1d
import h5py
import matplotlib.pyplot as plt


def traj_plotter(shotnumber, time, runnumber, detector_par, species, energy):
    rootFolder = os.path.abspath(os.getcwd())
    resultFolder = rootFolder+'/results'
    fieldFolder  = rootFolder+'/input/fieldGrid'
    cdbFolder    = rootFolder+'/input/cdb'
    ionProfileFolder =rootFolder+'/input/ionProf'

    shotAndTime = shotnumber + '_' + time

    workingFolder = resultFolder+'/'+shotAndTime+'/'+runnumber

    #ionisation
    ionR = numpy.genfromtxt(fname=ionProfileFolder+'/'+shotAndTime+'/'+species+'/'+energy+'/rad.dat')
    ionY = numpy.genfromtxt(fname=ionProfileFolder+'/'+shotAndTime+'/'+species+'/'+energy+'/ionyeald.dat')
    ionX1 = numpy.gradient(ionY)
    ionF = interp1d(ionR, -ionX1/numpy.amax(numpy.absolute(ionX1)))

    plt.plot(ionR, ionF(ionR))

    try:
        print('Save plot to '+workingFolder+'/ion_'+shotAndTime+'.pdf')
        plt.savefig(workingFolder+'/ion_'+shotAndTime+'.pdf')
        plt.clf()

    except:
        print('Unable to save to : '+workingFolder+'/ion_'+shotAndTime+'.pdf')
        plt.show()

    # trajectories
    t_rad = numpy.genfromtxt(fname=workingFolder+'/t_rad.dat')
    t_z   = numpy.genfromtxt(fname=workingFolder+'/t_z.dat')
    t_tor = numpy.genfromtxt(fname=workingFolder+'/t_tor.dat')
    cellid = numpy.genfromtxt(fname=workingFolder+'/detector/cellid.dat')

    rcord = numpy.genfromtxt(fname=fieldFolder+'/'+shotAndTime+'/rcord.dat')
    zcord = numpy.genfromtxt(fname=fieldFolder+'/'+shotAndTime+'/zcord.dat')
    psi2  = numpy.genfromtxt(fname=fieldFolder+'/'+shotAndTime+'/psi2.dat')

    hdf5Path = cdbFolder+'/'+shotnumber+'/EFITXX/EFITXX.1.h5'
    hdf5File = h5py.File(hdf5Path, 'r')
    hdf5Times = hdf5File.get('time')[()]
    hdf5PsiBoundary = hdf5File.get('output/globalParameters/psiBoundary')[()]
    hdf5PsiAxis = hdf5File.get('output/globalParameters/psiAxis')[()]
    hdf5BoundaryR = hdf5File.get('output/separatrixGeometry/boundaryCoordsR')[()]
    hdf5BoundaryZ = hdf5File.get('output/separatrixGeometry/boundaryCoordsZ')[()]

    hdf5Index = int(numpy.argwhere(hdf5Times == float(time)/1000.))

    psiAxis = hdf5PsiAxis[hdf5Index]
    psiBoundary = hdf5PsiBoundary[hdf5Index]
    boundaryR = hdf5BoundaryR[hdf5Index,:]
    boundaryZ = hdf5BoundaryZ[hdf5Index,:]

    PSI = numpy.transpose(psi2)/psiAxis

    R, Z =  numpy.meshgrid(rcord, zcord)

    t_index = cellid >= 0
    t_rad = t_rad[:,t_index]
    t_z   = t_z  [:,t_index]
    t_tor = t_tor[:,t_index]
    numberOfIons = numpy.size(t_rad, 1)

    #ion histogram (real)
    fig = plt.figure()
    plt.hist(t_rad[0,:], density=True, bins=50)
    plt.xlabel(r'$R$ [m]')
    plt.title (r'COMPASS #'+shotnumber+'\n ($t = $'+time+' ms)')

    try:
        print('Save plot to '+workingFolder+'/start_'+shotAndTime+'.pdf')
        plt.savefig(workingFolder+'/start_'+shotAndTime+'.pdf')
        plt.clf()

    except:
        print('Unable to save to : '+workingFolder+'/start_'+shotAndTime+'.pdf')
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

    fig.suptitle(r'COMPASS #'+shotnumber+' ($t = $'+time+' ms)')
    #plt.colorbar()

    try:
        print('Save plot to '+workingFolder+'/end_'+shotAndTime+'.pdf')
        plt.savefig(workingFolder+'/end_'+shotAndTime+'.pdf')
        plt.clf()

    except:
        print('Unable to save to : '+workingFolder+'/end_'+shotAndTime+'.pdf')
        plt.show()

    #2D figure
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #plt.contourf(R, Z, -PSI, levels=numpy.arange(-1.05,2,.05), cmap='Spectral')
    #plt.contourf(R, Z, PSI, levels=numpy.arange(-2,2,.05), cmap='coolwarm')
    #plt.contourf(R, Z, PSI, levels=numpy.arange(-2,2,.05), cmap='twilight_shifted')
    plt.contourf(R, Z, PSI, levels=numpy.arange(-2,2,.1), cmap='GnBu')
    plt.plot(boundaryR, boundaryZ, color=(0,0,0))
    for i in range(numberOfIons):
        i1 = 1
        i2 = numpy.asscalar(ionF(t_rad[1,i]))
        plt.plot([rcord[-1], t_rad[1,i]], [t_z[0,i], t_z[0,i]], color=(i1, i1, 0))
        if i2<0:
            i2=0
        if i2>1:
            i2=1
        plt.plot(t_rad[:,i], t_z[:,i], color=(i2, 0, 0))

    # plot detector
    try:
        detector_y_size = 12.2e-3
        detector_x0, detector_y0, detector_z0, detector_phi, detector_psi = map(float, detector_par.split(','))

        detector_x_tilt = detector_y_size * numpy.sin(detector_phi / 180 * numpy.pi)
        detector_y_tilt = detector_y_size * numpy.cos(detector_phi / 180 * numpy.pi)
        plt.plot([detector_x0, detector_x0], [detector_y0, max(zcord)], color=(0, 0, 1))
        plt.plot([detector_x0, detector_x0], [detector_y0, max(zcord)], color=(0, 0, 1))
        plt.plot([detector_x0 - detector_x_tilt, detector_x0 + detector_x_tilt],
                 [detector_y0 + detector_y_tilt, detector_y0 - detector_y_tilt],
                 color=(0, 0, 1), linewidth=2)
    except:
        print('Invalid detector coordinates.')
        pass

    plt.xlabel(r'$R$ [m]')
    plt.ylabel(r'$Z$ [m]')
    plt.title (r'COMPASS #'+shotnumber+'\n ($t = $'+time+' ms)')
    ax.set_aspect('equal', 'box')
    #plt.colorbar()

    try:
        print('Save plot to '+workingFolder+'/traj_'+shotAndTime+'.pdf')
        plt.savefig(workingFolder+'/traj_'+shotAndTime+'.pdf')
        plt.clf()

    except:
        print('Unable to save to : '+workingFolder+'/traj_'+shotAndTime+'.pdf')
        plt.show()

