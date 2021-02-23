"""
fig1_cluster_results.py

Plot the cluster results in Figure 1 of the paper.

created: Christopher Brittin
date: 22 February 2021

"""
import os
import sys
sys.path.append(r'./analysis')
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import matplotlib.pyplot as plt


from mc_community_plot_figure import *


CONFIG = 'configs/config.ini'
mpl.rcParams['ytick.labelsize'] = 6



def run(_cfg,fout=None,source_data=None):
    _doc = ("Three figures are produced. " 
            "Figure 1 is a cluster frequency plot."
            "Figure 2 is the cluster frequency plot with rows/cols reordered to match "
            "the A/P ordering. Figure 2 matches Fig. 1c in the paper. "
            "Figure 3 reproduce Fig. 1d in the paper.")
    print(_doc)
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)
    perturbations = 'data/perturbations/mc_cluster_rand_sig23_m4_t35.npz' 
    clusters = 'data/clusters/clusters_s23_m4_t35.csv'
    clusters = aux.read.into_dict(clusters)
    no_ticklabels = False
    reorder_nodes = True
    no_cbar = False
    #run(params.config,fout=params.fout,source_data = params.source_data)

    im = plot_clustermap(perturbations,cfg,clusters,
            no_cbar=no_cbar,no_ticklabels=no_ticklabels,reorder=reorder_nodes)

    plot_compare(cfg,im) 
    
    plt.show()
 
if __name__=="__main__":
    run(CONFIG) 
