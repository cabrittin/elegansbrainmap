"""
mht.py

Functions for doing multiple-hypothesis testing

Creator: Christopher Brittin
Date: 17 February 2016

"""

from scipy.stats import binom
from scipy.stats import norm
from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt


def pval_binom(data,p):
    """
    Computes p-values assuming binomial distribution.
    Input:
         data: Nx3 array
               [number of events, actual events, predicted events]
         p: probability defining distribution
    Output:
          pval: N array of pvals
    """
    
    pval = binom.cdf(data[:,1],data[:,0],p)
    idx = np.where(data[:,1] > data[:,2])
    pval[idx] = 1 - pval[idx]
    return pval

def pval_norm(data):
    """
    Computes p-values assuming standard normal distribution.
    Input:
         data: Nx3 array
               [test events, mu, variance]
    Output:
          pval: N array of pvals
    """
    x = (data[:,0] + 0.5 - data[:,1])/np.sqrt(data[:,2])
    pval = norm.cdf(x,loc=0,scale=1)
    idx = np.where(x > 0)
    pval[idx] = 1 - pval[idx]
    return pval    

def pval_norm_one_sided(data,**kwargs):
    """
    Computes p-values assuming standard normal distribution.
    Tests on
    Input:
         data: Nx3 array
               [test events, mu, variance]
         positive: (optional,boolean, default True) 
                   test for positive or negative side of mean
    Output:
          pval: N array of pvals
    """    

    positive = True
    if 'positive' in kwargs: positive=kwargs['side']
    x = (data[:,0] + 0.5 - data[:,1])/np.sqrt(data[:,2])
    pval = norm.cdf(x,loc=0,scale=1)
    if positive: pval = 1 - pval
    return pval

    

def adjusted_pval(pval):
    """
    Computes adjusted pvals using the 
    Benjamini-Hochberg procedure.
    
    Input:
         pval: N array of pvals

    Output:
         adjpval = float value
    """
    pval.sort()
    N = float(pval.shape[0])
    #for i in range(N):
    #    pval[i] = min(N*pval[i]/(i+1),1)
    idx = np.arange(N) + 1
    pval = N*pval/idx
    pval[pval > 1] == 1
    print(pval)
    adjpval = np.min(pval)
    return adjpval

def plot_actual_vs_predict(ax,data,colorbar=False,vmin=0,vmax=0.05,
                           tx=0.65,ty=0.1,th=0.05):
    """
    Plots actual vs predicted values
   
    Input:
        ax : handle for plot
        data: Nx2 array
              [neighbors,actual value, predicted value, variance]
    
    """    
    [N,m] = data.shape
    r =pearsonr(data[:,1],data[:,2])
    pval = pval_norm(data[:,1:])
    scat = ax.scatter(data[:,1],data[:,2],c=pval,s=72,
                      cmap=plt.get_cmap('jet'),
                      vmin=vmin,vmax=vmax)
    fit = np.max(data)
    fit = [0,fit]
    ax.plot(fit,fit,'r')
    ax.set_xlim(fit)
    ax.set_ylim(fit)
    ax.text(tx,ty+1.5*th,r'$n$ = %d'%N,
            transform=ax.transAxes,fontsize=24)
    ax.text(tx,ty,r'$r^2$ = %1.2f' %r[0]**2,
            transform=ax.transAxes,fontsize=24)
    if colorbar: cbar = plt.colorbar(scat)
    adj_pval = adjusted_pval(np.copy(pval))
    print(adj_pval)
    if adj_pval < 0.01:
        ax.text(tx,ty-th,r'$p_{adj}$ < 0.01',
                transform=ax.transAxes,fontsize=24)
    else:
        ax.text(tx,ty-th,r'$p_{adj}$=%1.2f' %adj_pval,
                transform=ax.transAxes,fontsize=24)

