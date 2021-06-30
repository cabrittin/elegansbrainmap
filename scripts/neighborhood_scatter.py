import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import defaultdict
import numpy as np
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                          mark_inset)

mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5

from connectome.load import from_db
import ioaux

REMOVE = ['VB01', 'VD01','VC01']
WIN_SIZE = 6
WIN_SIZE_HALF = int(WIN_SIZE / 2)
CONFIG = os.environ['CONFIG']
BAD = ['PLNL']

def run(_cfg,fout=None,source_data=None):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)

    left = ioaux.read.into_list(cfg['mat']['left_nodes'])
    right = ioaux.read.into_list(cfg['mat']['right_nodes'])
    cclass = ioaux.read.into_dict(cfg['mat']['class']) 

    N = from_db('N2U',adjacency=True,dataType='networkx',remove=REMOVE)
    J = from_db('JSH',adjacency=True,dataType='networkx',remove=REMOVE)
   
    SD = [['cell_class','adult_left','adult_right','l4_left','l4_right']]

    adult_left,adult_right,l4_left,l4_right = [],[],[],[]
    delta_l4,delta_adult = [],[]
    for (l,r) in zip(left,right):
        if l in BAD: continue
        if not N.A.has_node(l): continue
        if not N.A.has_node(r): continue
        if not J.A.has_node(l): continue
        if not J.A.has_node(r): continue
        adult_left.append(N.A.degree(l))
        adult_right.append(N.A.degree(r))
        l4_left.append(J.A.degree(l))
        l4_right.append(J.A.degree(r))
        delta_l4.append(abs(J.A.degree(l) - J.A.degree(r)))
        delta_adult.append(abs(N.A.degree(l) - N.A.degree(r)))
        tmp = [cclass[l],N.A.degree(l),N.A.degree(r),J.A.degree(l),J.A.degree(r)]
        SD.append(tmp)

    i = np.argsort(adult_left)
    delta_l4 = np.array(delta_l4)[i]
    delta_adult = np.array(delta_adult)[i]
    adult_left = np.array(adult_left)

    #Plot data
    fig,ax = plt.subplots(1,1,figsize=(3,3))
    
    s = 4
    print(len(adult_left))
    ax.scatter(adult_left,adult_left,color = 'k', marker='s',s=s,label="Adult left")
    ax.scatter(adult_left,adult_right,color='b',marker='o',s=s,label="Adult right")
    ax.scatter(adult_left,l4_left,color='g',marker='^',s=s,label="L4 left")
    ax.scatter(adult_left,l4_right,color='r',marker="v",s=s,label="L4 right")

    ax.set_xlabel('Adult left neighborhood size',fontsize=7)
    ax.set_ylabel('Neighborhood size',fontsize=7)
    ax.legend(fontsize=5)
    
    ax2 = plt.axes([0,0,1,1])
    ip = InsetPosition(ax, [0.67,0.13,0.3,0.3]) 
    ax2.set_axes_locator(ip)
    ax2.bar(adult_left[i],delta_l4,color='r',linewidth=3,label='L4')
    ax2.bar(adult_left[i],delta_adult,color='b',linewidth=3,label='Adult')
    ax2.set_xlabel('Size',fontsize=6)
    ax2.set_ylabel('Spread',fontsize=6)
    ax2.legend(fontsize=5)
    ax2.set_ylim([0,30])
    
    for tick in ax2.xaxis.get_major_ticks(): tick.label.set_fontsize(5) 
    for tick in ax2.yaxis.get_major_ticks(): tick.label.set_fontsize(5)
    plt.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()
    if source_data: ioaux.write.from_list(source_data,SD)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    
    parser.add_argument('-o','--fout',
                        action = 'store',
                        dest = 'fout',
                        default = None,
                        required = False,
                        help = 'Output svg file')
    
    parser.add_argument('-sd','--source_data',
                        action = 'store',
                        dest = 'source_data',
                        default = None,
                        required = False,
                        help = 'Output source data file')
    
    params = parser.parse_args()
    run(params.config,fout=params.fout,source_data = params.source_data)
