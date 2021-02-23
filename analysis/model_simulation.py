"""
@name: model_simulation.py
@description:
Script to generate K surrogate datasets for axonal, synaptic or
gap junction contacts, based on 3 parameter model
f = fraction of target contacts
p = probability to connect to target contact
q = 1-s = probability to connect to non-target contact.

Brittin et al 2020 
March 2020

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05
"""

import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from tabulate import tabulate

mpl.rcParams['ytick.labelsize'] = 5
mpl.rcParams['xtick.labelsize'] = 5

COLOR = '#a20000'

def get_params(cfg,contact):
    params = np.array([cfg.getfloat(contact,p) for p in ['f','p','s']])
    params[2] = 1 - params[2]
    return params

def run_simulation(K,n,params):
    n,K = int(n),int(K)
    nt = int(np.round(n*params[0]))
    no = n - nt

    C = np.zeros((n,K))
    Rt = np.random.rand(nt,K)
    Ro = np.random.rand(no,K)

    idx_t = np.where(Rt < params[1])
    idx_o = np.where(Ro < params[2])
    idx_o = (idx_o[0] + nt,idx_o[1])
    C[idx_t] = 1
    C[idx_o] = 1

    num_core = C[:nt,:].sum(axis=0)
    num_var = C[nt:,:].sum(axis=0)
    num_all = C.sum(axis=0)
    stats = np.array([[num_all.mean(),num_all.std()],
                      [num_core.mean(),num_core.std()],
                      [num_var.mean(),num_var.std()]])
    
    sumC = C.sum(axis=1)
    A = np.zeros(K+1)
    for k in range(K+1):
        A[k] = len(np.where(sumC==k)[0])

    return np.array([n,nt,no]),stats,A

def print_stats(N,stats):
    stats2 = stats / stats[0,0]
    stats = stats.tolist()
    stats2 = stats2.tolist()
    print(tabulate([N.tolist()],headers=['N','Nt','No'])) 
    print('\nCounts')
    label = ['f','ft','fo']
    tmp = [[f] + stats[i] for (i,f) in enumerate(label)]
    print(tabulate(tmp,headers=['mean','std']))
    print('\nNormalized')
    tmp = [[f] + stats2[i] for (i,f) in enumerate(label)]
    print(tabulate(tmp,headers=['mean','std']))

def plot_surrogate(ax1,A):
    n = A.shape[0]
    cumsum = np.cumsum(A)
    Anorm = float(cumsum[-1] - A[0])
    A[1:] /= Anorm
    ax2 = ax1.twinx()
    
    def rescale(x):
        return x * Anorm
    
    def unscale(x):
        return x / Anorm

    def scale_counts(ax1):
        y1,y2 = ax1.get_ylim()
        ax2.set_ylim(rescale(y1),rescale(y2))
        ax2.figure.canvas.draw()

    ax1.callbacks.connect("ylim_changed", scale_counts)
    ax1.plot(range(1,n),A[1:],color='k',linewidth=0.8)
    #ax1.set_xlim(0, n-0.5)
    xticks = [i for i in range(0,int(n),int(n/4))]
    #if n > 100: xticks = [i for i in range(0,int(n),int(n/2))]
    ax1.set_xticks(xticks)
    ax1.set_xlim(0,xticks[-1])
    ylim1 =  ax1.get_ylim()
    ax1.set_ylim(0,ylim1[1])
    yticks1 = [unscale(y) for y in ax2.get_yticks()]
    ax1.set_yticks(yticks1) 
    if max(yticks1) < 0.05:
        ytks = ['%1.3f'%y for y in yticks1]
    else:
        ytks = ['%1.2f'%y for y in yticks1]
    ax1.set_yticklabels(ytks)
    ylim2 = ax2.get_ylim()
    ax1.set_ylim(0,unscale(ylim2[1]))
    ax2.set_ylim(0,ylim2[1])
    ax1.spines['left'].set_color(COLOR)
    ax1.set_xlabel('$\delta$',fontsize=7)
    ax1.tick_params(axis='y', colors=COLOR)
    #ax1.set_ylabel('Fraction of contact sites ($\delta$)',fontsize=8)
    #ax2.set_ylabel('# of contact sites ($\delta$)',fontsize=8)
    return ax1,ax2

MODE = {'model_m':0,'model_c':1,'model_g':2}
fw,fh = 4.1,3.9
fw,fh = 3.5,3.4
fw,fh = 3.5,2.9
fw,fh = 3.2,2.9
CONFIG = './configs/config.ini'

#YTICKS1 = [[0.00,0.08,0.16,0.24,0.32,0.40],[0.00,0.03,0.06,0.09,0.12,0.15,0.18],[0.00,0.02,0.04,0.06,0.08],[0.000,0.007,0.014,0.021,0.028]]

YTICKS2 = [[0,250,500,750,1000],[0,100,200,300,400,500],[0,50,100,150,200,250],[0,20,40,60,80]]


YTICKS1 = {
        'model_m': [[0.00,0.08,0.16,0.24,0.32,0.40],[0.00,0.03,0.06,0.09,0.12,0.15,0.18],
                    [0.00,0.02,0.04,0.06,0.08],[0.000,0.007,0.014,0.021,0.028]],
        'model_c': [[0.00,0.08,0.16,0.24,0.32],[0.00,0.03,0.06,0.09,0.12,0.15],
                    [0.00,0.02,0.04,0.06,0.08],[0.000,0.007,0.014,0.021,0.028]],
        'model_g': [[0.00,0.07,0.14,0.21,0.28,0.35,0.42,],[0.00,0.03,0.06,0.09,0.12,0.15,0.18],
                    [0.00,0.02,0.04,0.06,0.08,0.10],[0.000,0.007,0.014,0.021,0.028,0.035,0.042]]
        }

def run(params):
    fout = params.fout
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
    
    m = params.mode
    D = np.load(cfg['model_data']['all'])
    fig,ax = plt.subplots(2,2,figsize=(fw,fh))
    _ax = ax.flatten()
    for (i,K) in enumerate([4,20,100,1000]):
        print('Number of surrogates: %d' %K)
        N = D.sum(axis=1)
        params = get_params(cfg,m)
        n,stats,A = run_simulation(K,N[MODE[m]],params)
        print_stats(n,stats)
        ax1,ax2 = plot_surrogate(_ax[i],A)
        ax1.set_yticks(YTICKS1[m][i])
        ax2.set_yticks(YTICKS2[i])
        tmp = '%1.2f'
        if i == 3: tmp='%1.3f'
        ystr = [tmp%_y for _y in YTICKS1[m][i]]
        ax1.set_yticklabels(ystr)
        if i in [0,2]: ax1.set_ylabel('Contact fraction',fontsize=7,color=COLOR)
        if i in [1,3]: ax2.set_ylabel('# of contacts',fontsize=7)
        #_ax[i].set_title('%d animals' %(K/2),fontsize=8)
        
    print(A)
    plt.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('mode',
                        action = 'store',
                        help = 'Mode: model_m, model_c, model_g')

    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    
    parser.add_argument('-o','--output',
                        dest = 'fout',
                        action = 'store',
                        default = None,
                        required = False,
                        help = 'Path to output file')
 
    params = parser.parse_args()
    run(params)

