import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as colormap
import scipy.stats as st


def detector(shotnumber, time, runnumber, is_debug=False):
    results_folder = os.path.join('results', f'{shotnumber}_{time}', runnumber)

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

        plt.xlabel(r"$x \mathrm{ [mm]}$")
        plt.ylabel(r"$y \mathrm{ [mm]}$")
        plt.title(f"COMPASS #{shotnumber} (t={time} s)")

        try:
            print(f'Save plot to {os.path.join(results_folder, f"detpy2_{shotnumber}_{time}.pdf")}')
            plt.savefig(os.path.join(results_folder, f"detpy2_{shotnumber}_{time}.pdf"))
            plt.clf()
        except Exception as e:
            print(f'Unable to save to: {results_folder} ({e})')
            plt.show()
    except Exception as e:
        print(f'Unable to plot: {results_folder} ({e})')
