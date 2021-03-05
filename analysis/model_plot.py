"""
@name: model_plot.py

@description
Plot the model fit

@author: Christopher Brittin
@email: 'cabrittin' + <at> + 'gmail' + '.' + 'com'
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['ytick.labelsize'] = 7
mpl.rcParams['xtick.labelsize'] = 7

from models.specificity import Model_P3_Avoidance as Avoidance

#THRESH = [0,0,0]
#VAL = [0.5,0.4,0.4]


def get_avoidance_error(error,val=0.5,high=False):
    A = Avoidance()
    aparams = A.min_error(error,derror=2000)
    params = None
    for i in range(2000):
        if high:
            if aparams[i,2] >= val:
                params = aparams[i,:3]
                break
        else:
            if aparams[i,2] <= val:
                params = aparams[i,:3]
                break
    return params

#def get_legend(params):
#    alabel = ('Model: ($p$,$s$,$f$) = (%1.2f,%1.2f,%1.2f)' 
#                %(params[0],1 - params[1],params[2]))
#    return alabel 

def get_legend(params):
    alabel = ('($p$,$s$,$f$) = (%1.2f,%1.2f,%1.2f)' 
                %(params[0],1 - params[1],params[2]))
    return alabel 


def plot_fit(ax,data,fit,ylabel=None,labels=None,legend=None,width=0.35,width2=0.33):
    n = data.sum()
    data /= n
    elabel = 'Empirical ($n$ = %d)' %n
    ind = np.arange(4)
    #ax.grid(zorder=1)
    ax.bar(ind-0.5*width,data,width=width2,color='k',zorder=2,label=elabel)
    ax.bar(ind+0.5*width,fit,width=width2,zorder=2,color='#a20000',label=legend)

    #ax.set_ylim([0,0.5])
    #x = [float('0.%d'%i) for i in range(6)]
    ax.set_ylim([0.0,0.8])
    x = [float('0.%d'%i) for i in range(0,9,2)]
    ax.set_yticks(x)
    ax.set_yticklabels(x,fontsize=5)
    ax.set_xlim([-0.5,3.5])
    ax.set_xticks(ind)
    if labels: ax.set_xticklabels(labels,fontsize=7)
    ax.legend(fontsize=6)
    #ax.set_ylabel('Fraction of contacts',fontsize=6) 
    ax.set_ylabel('Contact fraction',fontsize=7) 
    #if ylabel: ax.set_ylabel(ylabel,fontsize=6)

#figw,figh = 2.3,3.9
figw,figh = 1.9685,2.8

def run(params,fout=None,source_data=None):
    data = np.load(params.data)
    #data = data[:,1:]
    print(data)
    #data /= data.sum(axis=1)[:,None]
    avoidance = np.load(params.avoidance_error)
    print(data/data.sum(axis=1)[:,None])
    VAL = params.value
    THRESH = params.thresh

    A = Avoidance()
    _label = ['$\mathbb{%s}^1$','$\mathbb{%s}^2$','$\mathbb{%s}^3$','$\mathbb{%s}^4$']
    alabel = [l%'M' for l in _label]
    clabel = [l%'C' for l in _label]
    glabel = [l%'G' for l in _label]

    ddx = 0
    fig,ax = plt.subplots(3,1,figsize=(figw,figh))
    if params.fig_title: fig.canvas.set_window_title(params.fig_title)  
    try:
        print('Axon contacts')
        aparams = get_avoidance_error(avoidance['aerror'],val=VAL[0],high=THRESH[0])
        A.set(aparams[:3])
        leg = get_legend(aparams)
        print(leg)
        print(A.Z)
        plot_fit(ax[0],data[ddx,:],A.Z[1:],legend=leg,labels=alabel,ylabel='Fraction of axon contacts')
        #ax[0].set_title('Adjacency',fontsize=8)
        ax[0].legend(loc='upper left',fontsize=6)
        ddx += 1
    except:
        pass
    
    print('Chemical synapses')
    aparams = get_avoidance_error(avoidance['cerror'],val=VAL[1],high=THRESH[1])
    A.set(aparams[:3])
    leg = get_legend(aparams)
    print(leg)
    print(A.Z)
    plot_fit(ax[1],data[ddx,:],A.Z[1:],legend=leg, 
            labels=clabel,ylabel='Fraction of chemical syn.')
    ax[1].legend(loc='upper left',fontsize=6)
    #ax[1].set_title('Chemical syn.',fontsize=8) 
    ddx += 1

    print('Gap junctions')
    aparams = get_avoidance_error(avoidance['gerror'],val=VAL[2],high=THRESH[2])
    A.set(aparams[:3])
    leg = get_legend(aparams)
    print(leg)
    print(A.Z)
    plot_fit(ax[2],data[ddx,:],A.Z[1:],legend=leg,
            labels=glabel,ylabel='Fraction of gap junc.')
    #ax[2].set_title('Gap J.',fontsize=8)
    ax[2].legend(loc='upper left',fontsize=6)

    plt.tight_layout(pad=0.2)
    if params.fout: plt.savefig(params.fout)
    #plt.show()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('data',
                        action = 'store',
                        help = 'data file')
    
    parser.add_argument('avoidance_error',
                        action = 'store',
                        help = 'Avoidance error file')
    
    parser.add_argument('-o','--output',
                        action='store',
                        required=False,
                        dest='fout',
                        default=None,
                        help = 'Output file')
    
    params = parser.parse_args()

    run(params)
