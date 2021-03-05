"""
ed5_cluster_results.py

Plot the cluster results in ED Fig 5i-l of the paper.

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
    perturbations = 'data/perturbations/mc_cluster_rand_sig23_m4_t35.npz' 
    clusters = 'data/clusters/clusters_s23_m4_t35.csv'
    clusters = aux.read.into_dict(clusters)
    no_ticklabels = False
    reorder_nodes = False
    no_cbar = False
    
     
    pdir = 'data/perturbations/'
    pfile = ['mc_cluster_rand_sig23_m4_t35.npz','mc_cluster_rand_sig23_JSH.npz',
            'mc_cluster_rand_sig23_N2U.npz','mc_cluster_rand_sig0_m4.npz',
            'mc_cluster_rand_sig12_m4.npz','mc_cluster_rand_sig45_m4.npz',
            'mc_cluster_rand_sig90_m4.npz','mc_cluster_rand_sig23_m4_s1.npz',
            'mc_cluster_rand_sig23_m4_t0.npz']
    fig_title = ['ED5i M4', 'ED5i L4', 'ED5i Adult', 'ED5j sigma=0',
            'ED5j sigma=0.12','ED5j sigma=0.45', 'ED5j sigma=0.9',
            'ED5k Posterior lobed removed', 'ED5k Small contacts included']
    
    for (pf,ft) in zip(pfile,fig_title):
        perturbations = pdir + pf 
        im = plot_clustermap(perturbations,cfg,clusters,fig_title=ft,
                no_cbar=no_cbar,no_ticklabels=no_ticklabels,reorder=reorder_nodes)
        
    plt.show()
 
if __name__=="__main__":
    run(CONFIG) 
