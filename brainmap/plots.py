"""
@name: plots.py
@description:
Module to store plot functions

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-05
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import wilcoxon

def membrane_contacts(ax,Z,labels=None):
    cmap = sns.cm.rocket_r
    g = sns.heatmap(Z,ax=ax,vmax=1.,vmin=0,xticklabels=labels,yticklabels=labels,
            cbar_kws={'label': 'Fraction of membrane contact'},cmap=cmap)
    g.set_xticklabels(g.get_xticklabels(),rotation=30)
    #plt.tight_layout() 

def cell_bundle_contact(ax,z,bundles,bcolor):
    x = np.arange(len(bundles))
    ax.bar(x,z,width=0.5,color=bcolor)
    ax.set_xticks(x)
    ax.set_xticklabels(bundles,fontsize=8)
    ax.set_ylim([0,1])
    m = 1.0 / len(z)
    ax.axhline(m,color='k',linestyle='--',linewidth=2)
    ax.set_ylabel('Fraction of membrane contact',fontsize=12)
    #plt.tight_layout()

def cell_segment_contact(ax,Z,bundles,bcolor):
    [n,m] = Z.shape
    mu = np.mean(Z,axis=0)
    se = np.std(Z,axis=0) / np.sqrt(n)
    x = np.arange(len(bundles))
    ax.bar(x,mu,yerr=se,width=0.5,color=bcolor)
    ax.set_xticks(x)
    ax.set_xticklabels(bundles,fontsize=8)
    ax.set_ylim([0,1])
    m = 1.0 / len(mu)
    ax.axhline(m,color='k',linestyle='--',linewidth=2)
    ax.set_ylabel('Average segment contact',fontsize=12)
    

def cell_segment_contact_dist(ax,Z,bundles,bcolor,nbins=10):
    n = len(bundles)
    for (i,b) in enumerate(bundles):
        ax.hist(Z[:,i],bins=nbins,range=(0,1),histtype='step',color=bcolor[i],linewidth=2,label=b)

    ax.set_xlim([0,1])
    ax.legend(loc='upper center')
    ax.set_ylabel('# segments',fontsize=12)
    ax.set_xlabel('Fraction of mebrane contact',fontsize=12)

def intra_bundle_contact(ax,bid,Z,color='0.5',nbins=10):
    z = Z[:,bid] / Z.sum(axis=1)
    ax.hist(z,bins=nbins,range=(0,1),color=color,linewidth=1,edgecolor='k')
    ax.set_xlim([0,1])
    ax.set_ylabel('# segments', fontsize=12)
    ax.set_xlabel('Intrabundle contact',fontsize=12)

def boxplots(ax,data,labels=None,pval=[],showfliers=False,
                  annotscale=1.25,positions=None,width=0.6,
                  xlim=None,ylim=None,xlabel=None,ylabel=None,
                  title=None,yfs=32,xfs=32,tfs=32,pfs=24,fout=None,
                  boxwidth=2,capwidth=3,
                  colors = None):
    flierprops = dict(markersize=1,marker='d',markerfacecolor='k')
    medianprops = dict(linestyle='-',linewidth=0.5,color='k')
    whiskerprops = dict(linestyle='-',linewidth=0.3,color='k')
    boxprops = dict(linewidth=0.4,color='k',facecolor='#ABABAB')
    capprops = dict(linewidth=0.3)
    bp = ax.boxplot(data,positions=positions,vert=True,
                    patch_artist=True,
                    labels=labels,
                    medianprops=medianprops,
                    whiskerprops=whiskerprops,
                    boxprops=boxprops,
                    capprops=capprops,
                    showfliers=showfliers,
                    flierprops=flierprops,
                    widths=width)
    for i in range(len(data)):
        print(i,'(Lower,Upper)',(np.percentile(data[i],25),np.percentile(data[i],75)))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x',direction='out',labelsize=7)
    ax.tick_params(axis='y',length=0,labelsize=5)
    ax.grid(axis='y',color='0.9',linestyle='-',linewidth=1)
    if colors:
        for patch,color in zip(bp['boxes'],colors):
            patch.set_facecolor(color)        
    
    for (i,j,p) in pval:
        print(i,j,p,stars(p))
        cap1 = bp['caps'][2*i+1].get_ydata()
        cap2 = bp['caps'][2*j+1].get_ydata()
        y_max = np.max(np.concatenate((cap1,cap2)))
        y_min = np.min(np.concatenate((cap1,cap2)))
        ax.annotate("", xy=(i+1, y_max), xycoords='data',
                    xytext=(j+1, y_max), textcoords='data',
                    arrowprops=dict(arrowstyle="-", ec='k',
                                    connectionstyle="bar,fraction=0.2"))
        ax.text(0.5*(i+j+2), y_max+annotscale, stars(p),
                fontsize=pfs,
                horizontalalignment='center',
                verticalalignment='center')  

    if ylim: ax.set_ylim(ylim)
    if xlim: ax.set_xlim(xlim)
    if ylabel: ax.set_ylabel(ylabel,fontsize=yfs)
    ax.tick_params(axis='x',labelsize=xfs)
    if title: ax.set_title(title,fontsize=tfs,y=1.04)
    #if fout: plt.savefig(fout)
    return bp

def print_wilcoxon(data,label=None,alternative="two-sided"):
    stat,pval = wilcoxon(data,alternative=alternative)
    _tmp = (': n=%d,mu=%1.3f,std=%1.3f,se=%1.3f,p-val=%.2E'
           %(len(data),np.mean(data),
             np.std(data),np.std(data)/np.sqrt(len(data)),
             pval))
    if label:
        _tmp = label + _tmp
    print(_tmp)

def stars(p):
   if p < 0.0001:
       return "****"
   elif (p < 0.001):
       return "***"
   elif (p < 0.01):
       return "**"
   elif (p < 0.05):
       return "*"
   else:
       return "ns"

def plot_dist(ax,deg,fout=None,nbins=100,hrange=(0,1),
              xlim=None,ylim=None,density=False,cumulative=False,
              xlabel=None,ylabel=None,fs=38,linewidth=1):

    y,bins,_ = ax.hist(deg,bins=nbins,range=hrange,
                       histtype='step',density=density,
                       cumulative=cumulative,color='k',
                       linewidth=linewidth,label='empirical data')

    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    if xlabel: ax.set_xlabel(xlabel,fontsize=fs)
    if ylabel: ax.set_ylabel(ylabel,fontsize=fs)
    if fout: plt.savefig(fout)
      
