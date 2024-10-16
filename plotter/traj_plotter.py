import sys
import os
import numpy as np
import numpy.matlib
from scipy.interpolate import interp1d
import h5py
import matplotlib.pyplot as plt
import qrcode


def traj_plotter(shotnumber, time, runnumber, detector_par, species, energy, title=None, R_min=0.65):
    root_folder = os.path.abspath(os.getcwd())
    result_folder = os.path.join(root_folder, 'results')
    field_folder = os.path.join(root_folder, 'input', 'fieldGrid')
    cdb_folder = os.path.join(root_folder, 'input', 'cdb')
    ion_profile_folder = os.path.join(root_folder, 'input', 'ionProf')

    shot_and_time = f"{shotnumber}_{time}"
    working_folder = os.path.join(result_folder, shot_and_time, runnumber)

    if title is None:
        r'COMPASS #' + shotnumber + '\n ($t = $' + time + ' ms)'

    # Ionization
    ionR = np.genfromtxt(fname=os.path.join(ion_profile_folder, shot_and_time, species, energy, 'rad.dat'))
    ionY = np.genfromtxt(fname=os.path.join(ion_profile_folder, shot_and_time, species, energy, 'ionyeald.dat'))
    ionX1 = np.gradient(ionY)
    ionF = interp1d(ionR, -ionX1 / np.amax(np.absolute(ionX1[ionR > R_min])))

    plt.plot(ionR, ionF(ionR))

    try:
        print(f"Save plot to {os.path.join(working_folder, 'ion_' + shot_and_time + '.pdf')}")
        plt.savefig(os.path.join(working_folder, f"ion_{shot_and_time}.pdf"))
        plt.clf()
    except Exception as e:
        print(f"Unable to save to: {working_folder} ({e})")
        plt.show()

    # Trajectories
    t_rad = np.genfromtxt(fname=os.path.join(working_folder, 't_rad.dat'))
    t_z = np.genfromtxt(fname=os.path.join(working_folder, 't_z.dat'))
    t_tor = np.genfromtxt(fname=os.path.join(working_folder, 't_tor.dat'))
    cellid = np.genfromtxt(fname=os.path.join(working_folder, 'detector', 'cellid.dat'))

    rcord = np.genfromtxt(fname=os.path.join(field_folder, shot_and_time, 'rcord.dat'))
    zcord = np.genfromtxt(fname=os.path.join(field_folder, shot_and_time, 'zcord.dat'))
    psi2 = np.genfromtxt(fname=os.path.join(field_folder, shot_and_time, 'psi2.dat'))

    hdf5_path = os.path.join(cdb_folder, shotnumber, 'EFITXX', 'EFITXX.1.h5')
    hdf5_file = h5py.File(hdf5_path, 'r')
    hdf5_times = hdf5_file.get('time')[()]
    hdf5_psi_boundary = hdf5_file.get('output/globalParameters/psiBoundary')[()]
    hdf5_psi_axis = hdf5_file.get('output/globalParameters/psiAxis')[()]
    hdf5_boundary_r = hdf5_file.get('output/separatrixGeometry/boundaryCoordsR')[()]
    hdf5_boundary_z = hdf5_file.get('output/separatrixGeometry/boundaryCoordsZ')[()]
    hdf5_index = int(np.argwhere(hdf5_times == float(time) / 1000.))
    psi_axis = hdf5_psi_axis[hdf5_index]
    psi_boundary = hdf5_psi_boundary[hdf5_index]
    boundary_r = hdf5_boundary_r[hdf5_index, :]
    boundary_z = hdf5_boundary_z[hdf5_index, :]

    # Calculate PSI
    PSI = np.transpose(psi2) / psi_axis

    # Create meshgrid
    R, Z = np.meshgrid(rcord, zcord)

    # Filter data
    t_index = (cellid >= 0) & (t_rad[1, :] > R_min)
    t_rad = t_rad[:, t_index]
    t_z = t_z[:, t_index]
    t_tor = t_tor[:, t_index]
    number_of_ions = np.size(t_rad, 1)

    # Plot ion histogram
    fig, ax = plt.subplots()
    ax.hist(t_rad[0, :], density=True, bins=50)
    ax.set_xlabel(r'$R$ [m]')
    ax.set_title(title)

    # Save plot
    plot_filename = os.path.join(working_folder, f'start_{shot_and_time}.pdf')
    try:
        plt.savefig(plot_filename)
        plt.clf()
        print(f'Save plot to {plot_filename}')
    except:
        print(f'Unable to save to: {plot_filename}')
        plt.show()

    # Create detector plots
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(t_rad[-1, :], t_z[-1, :], '.')
    axs[0, 0].set_ylabel(r'$Z$ [m]')
    axs[0, 0].yaxis.set_ticks_position('both')

    axs[0, 1].plot(t_tor[-1, :], t_z[-1, :], '.')
    axs[0, 1].set_xlabel(r'$T$ [m]')
    axs[0, 1].set_ylabel(r'$Z$ [m]')
    axs[0, 1].yaxis.tick_right()
    axs[0, 1].yaxis.set_ticks_position('both')

    axs[1, 0].plot(t_rad[-1, :], t_tor[-1, :], '.')
    axs[1, 0].set_xlabel(r'$R$ [m]')
    axs[1, 0].set_ylabel(r'$T$ [m]')
    axs[1, 0].yaxis.set_ticks_position('both')

    axs[1, 1].axis('off')

    fig.suptitle(title)
    plt.show()

    # Create a 2D figure
    plt.clf()
    plt.rc('font', family='serif')
    plt.rc('text', usetex=True)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot filled contours
    plt.contourf(R, Z, PSI, levels=np.arange(-2, 2, 0.1), cmap='GnBu')
    #plt.contourf(R, Z, PSI, levels=np.arange(0, 12, 0.3), cmap='GnBu')

    # Plot boundary line
    plt.plot(boundary_r, boundary_z, color='black')

    # Iterate over ions
    intensity_indices = np.argsort(ionF(t_rad[1, :]))
    ion_factor = ionF(t_rad[1, :])
    for i in intensity_indices:
        # Adjust ion factor within bounds
        ion_factor[i] = max(0, min(1, ion_factor[i]))
        # Plot ion trajectory
        plt.plot(t_rad[:, i], t_z[:, i], color=(ion_factor[i], 0, 0))

    alpha_factor = 50/number_of_ions
    for i in intensity_indices:
        # Plot atom trajectory
        plt.plot([rcord[-1], t_rad[1, i]], [t_z[0, i], t_z[0, i]], color=(1, 1, 0), alpha=ion_factor[i] * alpha_factor)
    try:
        # Detector parameters
        detector_y_size = 12.2e-3
        detector_x0, detector_y0, detector_z0, detector_phi, detector_psi = map(float, detector_par.split(','))

        # Calculate tilted coordinates
        detector_x_tilt = detector_y_size * np.sin(detector_phi / 180 * np.pi)
        detector_y_tilt = detector_y_size * np.cos(detector_phi / 180 * np.pi)

        # Plot detector lines
        detector_colour = (0, 0.5, 0)
        plt.plot([detector_x0, detector_x0], [detector_y0, max(zcord)],
                 color=detector_colour, linewidth=1)
        plt.plot([detector_x0 - detector_x_tilt, detector_x0 + detector_x_tilt],
                 [detector_y0 + detector_y_tilt, detector_y0 - detector_y_tilt],
                 color=detector_colour, linewidth=2)
    except ValueError:
        print('Invalid detector coordinates.')

    plt.xlabel(r'$R$ [m]')
    plt.ylabel(r'$Z$ [m]')
    plt.title(title)
    ax.set_aspect('equal', 'box')
    qr_content = f'https://taiga-project.github.io?{runnumber}'
    print(f'Add QR code: {qr_content}')
    qr_image = qrcode.make(qr_content)
    ax_qr = fig.add_axes([0.24, 0.495, 0.08, 0.5], anchor='NE', zorder=10)
    ax_qr.imshow(qr_image, cmap='gray')
    ax_qr.axis('off')
    try:
        plot_filename = os.path.join(working_folder, f'traj_{runnumber}')
        plt.savefig(plot_filename+'.pdf', dpi=300, bbox_inches='tight')
        print(f'Save plot to {plot_filename}.pdf')
        plt.savefig(plot_filename+'.png', dpi=300, bbox_inches='tight')
        print(f'Save plot to {plot_filename}.png')
        plt.clf()
    except:
        print(f'Unable to save to: {plot_filename}.pdf/png')
        plt.show()
