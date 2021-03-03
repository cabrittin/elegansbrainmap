"""
ed4_model_simulation.py

Model simulation in ED 4de of the paper.

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
    print("Chemical synapses simulation")
    Params = namedtuple('params',['config','mode','fout'])
    params = Params(_cfg,'model_c',None)
    model_simulation.run(params)
    
    print("Gap junction simulation")
    Params = namedtuple('params',['config','mode','fout'])
    params = Params(_cfg,'model_g',None)
    model_simulation.run(params)
 
if __name__=="__main__":
    run(CONFIG) 
