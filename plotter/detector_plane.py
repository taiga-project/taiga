import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import qrcode


def detector_plane(shotnumber, time, runnumber, detector_par, title=None):
    max_distance = 1e-4
    detector_par_list = detector_par.split(',')

    try:
        R0 = float(detector_par_list[0])
    except ValueError:
        R0 = 0.685

    try:
        Z0 = float(detector_par_list[1])
    except ValueError:
        Z0 = 0.23

    try:
        T0 = float(detector_par_list[2])
    except ValueError:
        T0 = 0.0

    try:
        detector_RZ_angle = float(detector_par_list[3])
    except ValueError:
        detector_RZ_angle = 38.0

    try:
        detector_ZT_angle = float(detector_par_list[4])
    except ValueError:
        detector_ZT_angle = 0.0

    results_folder = os.path.join('results', f'{shotnumber}_{time}', runnumber)

    if title is None:
        r'COMPASS #' + shotnumber + '\n ($t = $' + time + ' ms)'

    try:
        R = np.loadtxt(os.path.join(results_folder, 'rad.dat'))
        Z = np.loadtxt(os.path.join(results_folder, 'z.dat'))
        T = np.loadtxt(os.path.join(results_folder, 'tor.dat'))
    except FileNotFoundError:
        print(f'Invalid input folder: {results_folder}')
        return

    try:
        particle_on_detector = np.abs((Z - Z0) * np.tan((detector_RZ_angle / 180.0) * np.pi) + (R - R0)) < max_distance
        x = T[particle_on_detector]
        y = Z[particle_on_detector]
        xmin, xmax = -4e-2, 4e-2
        ymin, ymax = 18e-2, 26e-2
        #ymin, ymax = 24e-2, 32e-2
        #xmin, xmax = -8e-2, 0e-2
        #ymin, ymax = 14e-2, 22e-2

        # Perform the kernel density estimate
        xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([x, y])
        kernel = st.gaussian_kde(values, 0.05)
        f = np.reshape(kernel(positions).T, xx.shape)

        fig, ax = plt.subplots()
        plt.rc('font', family='serif')
        plt.rc('text', usetex=True)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect('equal')

        # Create a 2D histogram plot
        ax.hist2d(x, y, bins=250, range=[[xmin, xmax], [ymin, ymax]], cmap='afmhot')

        # Set labels
        plt.xlabel(r"$T$ [m]")
        plt.ylabel(r"$Z$ [m]")
        plt.title(title)

        qr_content = f'raw.githubusercontent.com/taiga-project/refs/run/{runnumber}'
        print(f'Add QR code: {qr_content}')
        qr_image = qrcode.make(qr_content)
        ax_qr = fig.add_axes([0.24, 0.495, 0.08, 0.5], anchor='NE', zorder=10)
        ax_qr.imshow(qr_image, cmap='gray')
        ax_qr.axis('off')

        try:
            print(f'Save plot to {os.path.join(results_folder, f"detpy_{shotnumber}_{time}.pdf")}')
            plt.savefig(os.path.join(results_folder, f"detpy_{shotnumber}_{time}.pdf"), dpi=300, bbox_inches='tight')
            plt.clf()
        except Exception as e:
            print(f'Unable to save to: {results_folder} ({e})')
            plt.show()
    except Exception as e:
        print(f'Unable to plot: {results_folder} ({e})')
