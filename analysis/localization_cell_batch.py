"""
@name: localization_cell_batch.py
@description:
Compute spatial reproduciblity of contacts.
Used for ED 1d and Supplementary Information 2


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
from tqdm import tqdm

import aux
from localization_prob import *

#CONFIG = os.environ['CONFIG']
CONFIG = './configs/config.ini'


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
 
    params = parser.parse_args()
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
    left = aux.read.into_list(cfg['mat']['left_nodes'])
    right = aux.read.into_list(cfg['mat']['right_nodes'])
    lrmap = aux.read.into_lr_dict(cfg['mat']['lrmap'])
    remove = aux.read.into_list(cfg['mat']['remove'])
    cclass = aux.read.into_dict(cfg['mat']['class'])
    
    A4 = nx.read_graphml(cfg['refgraphs']['adj']%4)
    data = aux.read.into_list2(cfg['adj_align']['fout'])
    
    edict = {}
    for (i,(a,b)) in enumerate(A4.edges()):
        edict[i] = (a,b)
        A4[a][b]['idx'] = i

    left = [n for n in left if A4.has_node(n)]
    
    for cell in tqdm(left,'Cell: '):
        D = build_loc_data(cell,data,A4) 
        neigh = sorted([v for v in A4.neighbors(cell)])
        f = prob_response(D)
        fig,ax = plt.subplots(1,2,figsize=(8,2.5))
        plot_order2(ax[0],D)
        ax[0].set_title(f'Localization of {cell} process')
        plot_contact_rate2(ax[1],f)
        plt.tight_layout()
        plt.savefig('results/adj_loc/%s_localization.svg'%cell)
        plt.close()
