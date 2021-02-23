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

import model_simulation

CONFIG = 'configs/config.ini'

def run(_cfg,fout=None,source_data=None):
    Params = namedtuple('params',['config','mode','fout'])
    params = Params(_cfg,'model_m',None)
    model_simulation.run(params)
 
if __name__=="__main__":
    run(CONFIG) 
