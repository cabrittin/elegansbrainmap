"""
@name: model_fit.py

@description
Fit model data


@author: Christopher Brittin
@email: 'cabrittin' + <at> + 'gmail' + '.' + 'com'
"""


import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation   
import numpy as np
from tqdm import tqdm

from models.specificity import *

CONFIG = os.environ['CONFIG']

def run_model(D,P,order=3):
    n = order + 1
    m = 101**order
    i = 0 
    aerror = np.zeros([m,n]) 
    for p in tqdm(P.loop(D[0,:]),desc = 'A: error'):
        aerror[i,:] = p
        i += 1
    
    i = 0
    cerror = np.zeros([m,n]) 
    for p in tqdm(P.loop(D[1,:]),desc='C: error'):
        cerror[i,:] = p
        i += 1
    
    i = 0 
    gerror = np.zeros([m,n]) 
    for p in tqdm(P.loop(D[2,:]),desc='G: error'):
        gerror[i,:] = p
        i += 1

    return [aerror,cerror,gerror]

def run_avoidance(D,fout):
    print('Avoidance model')
    error = run_model(D,Model_P3_Avoidance())
    eout = fout.replace('.npz','_fit.npz')
    np.savez(eout,aerror=error[0],cerror=error[1],gerror=error[2])

def load_data(data):
    D = np.load(data)
    D /= D.sum(axis=1)[:,None]
    D2 = np.zeros((3,5))
    D2[:,1:] = D
    D = D2
    return D

def run(params):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
   
    print('Weight fit')
    ldin = cfg['model_data']['low_weight']
    hdin = cfg['model_data']['high_weight']
    din = [ldin,hdin]
    for d in din:
        print('processing: %s' %d)
        D = load_data(d)
        print(D)
        fout = d.replace('.npy','.npz')
        run_avoidance(D,fout)
    
    print('Intra-/Inter bundle fit')
    ldin = cfg['model_data']['inter_bundle']
    hdin = cfg['model_data']['intra_bundle']
    din = [ldin,hdin]
    for d in din:
        print('processing: %s' %d)
        D = load_data(d)
        print(D)
        fout = d.replace('.npy','.npz')
        run_avoidance(D,fout)
    print('All data fit')
    d = cfg['model_data']['all']
    print('processing: %s' %d)
    D = load_data(d)
    print(D)
    fout = d.replace('.npy','.npz')
    run_avoidance(D,fout)

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

    run(parmams)
