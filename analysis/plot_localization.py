"""
@name: plot_localization.py
@description:
Plot localization results

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05
"""

import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.patches as mpatches
from matplotlib.patches import PathPatch
import pandas as pd
import seaborn as sns

import ioaux
from localization import *

CONFIG = os.environ['CONFIG']
width = 3
height = 2

fin = ['data/localization/localization_nbins11.npz',
        'data/localization/localization_nbins26.npz',
        'data/localization/localization_nbins51.npz']
label = [3.6,1.4,0.7]
binsize = ['3.6 μm', '1.4 μm', '0.7 μm']

fout = 'results/loc_summary/localization_%s.svg'


def adjust_box_widths(g, fac):
    """
    Adjust the withs of a seaborn-generated boxplot.
    """

    # iterating through Axes instances
    for ax in g.axes:

        # iterating through axes artists:
        for c in ax.get_children():

            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5*(xmin+xmax)
                xhalf = 0.5*(xmax - xmin)

                # setting new width of box
                xmin_new = xmid-fac*xhalf
                xmax_new = xmid+fac*xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new

                # setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])

def plot_data(fig,ax,df,bp_y):
    col = ['r','#ffb311','b', 'g']
    hatch = ['','xxxx','----']
    chatch = []
    for c in col:
        for h in hatch:
            chatch.append((c,h))
    circ1 = mpatches.Patch( facecolor='0.5',label='bin = 3.6 $\mu$m')
    circ2= mpatches.Patch( facecolor='0.5',hatch=r'xxxx',label='bin = 1.4 $\mu$m',edgecolor='w')
    circ3 = mpatches.Patch(facecolor='0.5',hatch=r'||||',label='bin = 0.7 $\mu$m',edgecolor='w')
    flierprops = dict(markersize=1,marker='d',markerfacecolor='k')
    medianprops = dict(linestyle='-',linewidth=0.5,color='k')
    whiskerprops = dict(linestyle='-',linewidth=0.3,color='k')
    capprops = dict(linewidth=0.3)
    bp = sns.boxplot(x='#_datasets',y=bp_y,hue='bin_size',data=df,ax=ax,
                width=0.8,linewidth=0.3,
                flierprops=flierprops,medianprops=medianprops,capprops=capprops)
    for i, patch in enumerate(bp.artists):
        # Boxes from left to right
        #hatch = next(hatches)
        patch.set_hatch(chatch[i][1])
        patch.set_facecolor(chatch[i][0])
        patch.set_edgecolor('w')
    
    bp.legend(handles=[circ1,circ2,circ3],loc='upper left',fontsize=6)
    #ax.legend(fontsize=6)
    ax.set_ylim([0,1])
    #ax.set_ylabel('Fraction $\hat{z}$ contacts',fontsize=7)
    ax.set_xticks([0,1,2,3])
    ax.set_xticklabels([1,2,3,4],fontsize=6)
    for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(6) 
    #ax.set_xlabel('Number of datasets',fontsize=7)
    adjust_box_widths(fig, 0.8)
    plt.tight_layout() 

def run(**args):
    data = []
    fig,ax = plt.subplots(1,1,figsize=(width,height))
    for (i,_fin) in enumerate(fin):
        cells = np.load(_fin)['cells']
        _f = np.load(_fin)['freq']
        f = _f / _f.sum(axis=1,keepdims=True)
        ncols = f.shape[1]
        for (j,c) in enumerate(cells):
            for _k in range(ncols):
                k = _k + 1
                data.append([c,binsize[i],k,_f[j,_k],f[j,_k]])
    df = pd.DataFrame(data=data,columns=['cell','bin_size','#_datasets',
        'raw_count_z_contacts','frac_z_contacts'])
    plot_data(fig,ax,df,'frac_z_contacts')
    ax.set_ylabel('Fraction $\hat{z}$ contacts',fontsize=7)
    ax.set_xlabel('Number of datasets',fontsize=7)
    #plt.savefig(fout%'_frac_z')
    #df.to_csv('source_data/source_data_extended_data_1e.csv',index=False)
    
    data = []
    fig,ax = plt.subplots(1,1,figsize=(width,height))
    for (i,_fin) in enumerate(fin):
        cells = np.load(_fin)['cells']
        _f = np.load(_fin)['counts']
        f = _f / _f.sum(axis=1,keepdims=True)
        ncols = f.shape[1]
        for (j,c) in enumerate(cells):
            for _k in range(ncols):
                k = _k + 1
                data.append([c,binsize[i],k,_f[j,_k],f[j,_k]])
    df = pd.DataFrame(data=data,columns=['cell','bin_size','#_datasets',
        'raw_count_max_z','frac_max_z'])
    plot_data(fig,ax,df,'frac_max_z')
    ax.set_ylabel('Fraction $\max(\delta)_\hat{z}$ contacts',fontsize=7)
    ax.set_xlabel('Number of datasets',fontsize=7)
    #plt.savefig(fout%'_counts_z')
    #df.to_csv('source_data/source_data_extended_data_1f.csv',index=False)
 
    data = []
    fig,ax = plt.subplots(1,1,figsize=(width,height))
    for (i,_fin) in enumerate(fin):
        cells = np.load(_fin)['cells']
        _f = np.load(_fin)['cs4']
        fsum = _f.sum(axis=1, keepdims=True)
        fsum[fsum==0] = 1
        f = _f / fsum
        ncols = f.shape[1]
        for (j,c) in enumerate(cells):
            for _k in range(ncols):
                k = _k + 1
                data.append([c,binsize[i],k,_f[j,_k],f[j,_k]])
    df = pd.DataFrame(data=data,columns=['cell','bin_size','#_datasets',
        'raw_count_max_z_c4','frac_max_z_c4'])
    plot_data(fig,ax,df,'frac_max_z_c4')
    ax.set_ylabel('Fraction $\max(\delta)_\hat{z}$ contacts',fontsize=7)
    ax.set_xlabel('Number of datasets',fontsize=7)
    #plt.savefig(fout%'_c4_zmaz')
    #df.to_csv('source_data/source_data_extended_data_1h.csv',index=False)
    
    plt.show()


if __name__=="__main__":
    run()
