"""
fig3_cluster_results.py

Generate Figure 3a-c of the paper.

created: Christopher Brittin
date: 22 February 2021

"""
import os
import sys
sys.path.append(r'./analysis')
from configparser import ConfigParser,ExtendedInterpolation
import argparse
from collections import namedtuple

import c4_contact_matrix


CONFIG = 'configs/config.ini'


def run(_cfg,fout=None,source_data=None):
    perturbations = 'data/perturbations/mc_cluster_rand_sig23_m4_t35.npz' 
    clusters = 'data/clusters/clusters_s23_m4_t35.csv'
    no_ticklabels = False
    no_cbar = False
 
    Params = namedtuple('params',['config','fin','clusters','no_cbar','no_ticklabels'])
    params = Params(_cfg,perturbations,clusters,no_cbar,no_ticklabels)
    
    c4_contact_matrix.run(params)
 
if __name__=="__main__":
    run(CONFIG) 
