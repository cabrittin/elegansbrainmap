"""
ed3_model_results.py

Plot modelling results for ED Fig 3.

created: Christopher Brittin
date: March 2021

"""
import os
import sys
sys.path.append(r'./analysis')
from configparser import ConfigParser,ExtendedInterpolation
from collections import namedtuple
import argparse
import matplotlib.pyplot as plt

import model_plot

CONFIG = 'configs/config.ini'

def run(_cfg,fout=None,source_data=None):
    mdir = './data/model/'
    Params = namedtuple('params',['config','data','avoidance_error','fout','fig_title','thresh','value'])
    data_files = ['model_low_weight.npy','model_high_weight.npy',
            'model_inter_bundle.npy','model_intra_bundle.npy',
            'model_complete.npy','model_rm_1EM.npy','model_prune_polyad.npy',
            'model_white.npy','model_witvliet.npy']
    fit_files = ['model_low_weight_fit.npz','model_high_weight_fit.npz',
            'model_inter_bundle_fit.npz','model_intra_bundle_fit.npz',
            'model_complete_fit.npz','model_rm_1EM_fit.npz','model_prune_polyad_fit.npz',
            'model_white_fit.npz','model_witvliet_fit.npz']
    fig_titles = ['Below average contact area','Above average contact area',
            'Inter-cluster contacts','Intra-cluster contacts',
            'Complete dataset','1 EM section removed','Variable polyads removed',
            'White et al. 1986','Witvliet et al. 2020']
    thresh = [[0,0,0],[1,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,1],[0,0,0]]
    values = [[0.5,0.4,0.4],[0.4,0.4,0.4],[0.5,0.4,0.4],[0.5,0.4,0.4],[0.5,0.4,0.4],
            [0.5,0.4,0.4],[0.5,0.4,0.4],[0.4,0.4,0.4],[0.5,0.4,0.4]]

    for (df,ff,ft,th,val) in zip(data_files, fit_files, fig_titles,thresh,values):
        params = Params(_cfg,mdir + df, mdir + ff, None, ft, th, val)
        print('\nModel: ' + params.fig_title)
        model_plot.run(params)
    plt.show() 
if __name__=="__main__":
    run(CONFIG) 
