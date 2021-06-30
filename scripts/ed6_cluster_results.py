"""
ed6_cluster_results.py

Plot the cluster results in ED Fig 6 of the paper.

created: Christopher Brittin
date: 22 February 2021

"""
import os
import sys
sys.path.append(r'./analysis')
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import matplotlib.pyplot as plt


from cluster_population_plot_figure import *


CONFIG = 'configs/config.ini'
mpl.rcParams['ytick.labelsize'] = 6



def run(_cfg,fout=None,source_data=None):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)
    clusters = 'data/clusters/clusters_s23_m4_t35.csv'
    clusters = ioaux.read.into_dict(clusters)
    no_ticklabels = False
    reorder_nodes = False
    no_cbar = False
    
     
    pdir = 'data/perturbations/'
    pfile = ['mc_cluster_rand_sig0_m4.npz','mc_cluster_rand_sig0_m3.npz',
            'mc_cluster_rand_sig0_m2.npz','mc_cluster_rand_sig0_m1.npz',
            'mc_cluster_rand_sig0_JSH.npz','mc_cluster_rand_sig0_N2U.npz']
    
    fig_title = ['M4','M3','M2','M1','L4','Adult']
    
    for (pf,ft) in zip(pfile,fig_title):
        perturbations = pdir + pf 
        im = plot_clustermap(perturbations,cfg,clusters,fig_title=ft,
                no_cbar=no_cbar,no_ticklabels=no_ticklabels,reorder=reorder_nodes)
        
    plt.show()
 
if __name__=="__main__":
    run(CONFIG) 
