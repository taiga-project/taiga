import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as colormap
import scipy.stats as st
import qrcode


def detector(shotnumber, time, runnumber, title=None, is_debug=False):
    results_folder = os.path.join('results', f'{shotnumber}_{time}', runnumber)

    if title is None:
        r'COMPASS #' + shotnumber + '\n ($t = $' + time + ' ms)'

    try:
        cell_counts = np.loadtxt(os.path.join(results_folder, 'detector', 'cellcounter.dat'))
        det_x = np.loadtxt(os.path.join(results_folder, 'detector', 'detx'))
        det_y = np.loadtxt(os.path.join(results_folder, 'detector', 'dety'))
    except FileNotFoundError:
        print(f'Invalid input folder: {results_folder}')
        return

    xmin, xmax = det_x[0, 0], det_x[-1, -1]
    ymin, ymax = det_y[0, 0], det_y[-1, -1]
    xmean = (det_x[:, 0] + det_x[:, 1]) / 2
    ymean = (det_y[:, 0] + det_y[:, 1]) / 2
    cmax = cell_counts.max() if cell_counts.max() > 0 else 1

    try:
        plt.rc('font', family='serif')
        plt.rc('text', usetex=True)
        fig, ax = plt.subplots()
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect('equal')

        # Contourf plot
        plt.contourf(xmean, ymean, cell_counts, cmap='afmhot', corner_mask=False)

        # Detector cell matrix plot
        for i in range(det_x[:, 0].size):
            for j in range(det_y[:, 0].size):
                rect = patches.Rectangle((det_x[i, 0], det_y[j, 0]),
                                         det_x[i, 1] - det_x[i, 0],
                                         det_y[j, 1] - det_y[j, 0],
                                         linewidth=1,
                                         edgecolor='k',
                                         facecolor=colormap.afmhot(cell_counts[j, i] / cmax))
                ax.add_patch(rect)

        if is_debug:
            print(cell_counts)

        plt.xlabel(r"toroidal detector position [mm]")
        plt.ylabel(r"vertical detector position [mm]")
        plt.title(title)

        qr_content = f'raw.githubusercontent.com/taiga-project/refs/run/{runnumber}'
        print(f'Add QR code: {qr_content}')
        qr_image = qrcode.make(qr_content)
        ax_qr = fig.add_axes([0.175, 0.495, 0.08, 0.5], anchor='NE', zorder=10)
        ax_qr.imshow(qr_image, cmap='gray')
        ax_qr.axis('off')

        try:
            filename = f"detector_{runnumber}.pdf"
            print(f'Save plot to {os.path.join(results_folder, filename)}')
            plt.savefig(os.path.join(results_folder, filename), dpi=300, bbox_inches='tight')
            plt.clf()
        except Exception as e:
            print(f'Unable to save to: {results_folder} ({e})')
            plt.show()
    except Exception as e:
        print(f'Unable to plot: {results_folder} ({e})')
