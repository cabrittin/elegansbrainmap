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


import model_plot


CONFIG = 'configs/config.ini'

def run(_cfg,fout=None,source_data=None):
    Params = namedtuple('params',['data','avoidance_error','fout'])
    params = Params('./data/model/model.npy',
                    './data/model/model_fit.npz',
                    None)
    
    model_plot.run(params)
 
if __name__=="__main__":
    run(CONFIG) 
