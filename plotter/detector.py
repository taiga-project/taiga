import sys
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.patches as patches
import matplotlib.cm as colormap
import scipy.stats as st


def detector(shotnumber, time, runnumber, is_debug=False):

    results_folder = 'results/'+shotnumber+'_'+time+'/'+runnumber+'/'
    
    try:
        C = np.loadtxt(results_folder+'detector/cellcounter.dat')
        X = np.loadtxt(results_folder+'detector/detx')
        Y = np.loadtxt(results_folder+'detector/dety')
    except:
        print ('Invalid input folder: '+results_folder)
    
    
    xmin = X[0,0]
    xmax = X[-1,-1]
    ymin = Y[0,0]
    ymax = Y[-1,-1]
    xmean = (X[:,0]+X[:,1])/2
    ymean = (Y[:,0]+Y[:,1])/2
    cmax = C.max()

    if cmax == 0:
        cmax = 1
    
    try:
        fig = pl.figure()
        ax = fig.gca()
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect('equal')
    
        # Contourf plot
        pl.contourf(xmean, ymean, C, cmap='afmhot', corner_mask=False)
    
        # Detector cell matrix plot
        for i in range (X[:,0].size):
            for j in range (Y[:,0].size):
                rect = patches.Rectangle((X[i,0], Y[j,0]), X[i,1]-X[i,0], Y[j,1]-Y[j,0], linewidth=1,
                                         edgecolor='k', facecolor=colormap.afmhot(C[j,i]/cmax))
                ax.add_patch(rect)

        if is_debug:
            print(C)
    
        pl.xlabel(r"$x \mathrm{ [mm]}$")
        pl.ylabel(r"$y \mathrm{ [mm]}$")
        pl.title(r"COMPASS #"+shotnumber+" (t="+time+" s)")
    
        try:
            print('Save plot to '+results_folder+'detpy2_'+shotnumber+'_'+time+'.pdf')
            pl.savefig(results_folder+'detpy2_'+shotnumber+'_'+time+'.pdf')
            pl.clf()
        except:
            print('Unable to save to : '+results_folder)
            pl.show()
    except:
        print('Unable to plot: '+results_folder)

