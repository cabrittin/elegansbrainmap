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
from collections import namedtuple
import argparse
import matplotlib.pyplot as plt

import model_init
import model_plot
import model_fit

CONFIG = 'configs/config.ini'

def run(_cfg,fout=None,source_data=None):
    mdir = './data/model/'
    isDir = os.path.isdir(mdir)
    if not isDir:
        print("%s directory does not exist....creating"%mdir)
        os.mkdir(mdir)
    Params = namedtuple('params',['config','delta','data','avoidance_error','fout',
                                    'fig_title','thresh','value'])
    params = Params(_cfg,4,
                    mdir + 'model.npy',
                    mdir + 'model_fit.npz',
                    None,None,[0,0,0],[0.5,0.4,0.4])
    
    isFile = os.path.isfile(params.data)
    if not isFile:
        print("Model data (%s) not found....."%params.data)
        print("Generating model data, this will only need to be done once.")
        model_init.run(params)
    
    isFile = os.path.isfile(params.avoidance_error)
    if not isFile:
        print("Model fit error (%s) not found....."%params.avoidance_error)
        print("Fitting model for all parameter combinations.")
        print("This will take some time, but this will only need to be done once.")
        model_fit.run(params)
    
    model_plot.run(params)
    plt.show()
 
if __name__=="__main__":
    run(CONFIG) 
