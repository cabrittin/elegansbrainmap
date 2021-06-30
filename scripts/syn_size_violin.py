"""
@name: syn_size_violoin.py
@description:
    plot distribution of synapse sizes

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-03
"""

import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator

from connectome.load import from_db
from connectome.format_graphs import *
#from networks.classification_measures import *
import ioaux

mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'

def clean_graph(Ref,A,idx=4):
    H = nx.Graph()
    if Ref.is_directed(): H = nx.DiGraph()
    for (a,b) in Ref.edges():
        if A.has_edge(a,b) and Ref[a][b]['id'] <= idx:
            H.add_edge(a,b,weight=Ref[a][b]['weight'],id=Ref[a][b]['id'])
    return H

def build_data(Ref,label):
    data = []
    for (a,b) in Ref.edges():
        data.append([Ref[a][b]['id'],Ref[a][b]['weight'],label])
    return data

def cohend(x,y):
    mu1 = np.mean(x)
    mu2 = np.mean(y)
    n1 = len(x) - 1
    n2 = len(y) - 1
    sig1 = np.std(x,ddof=1)
    sig2 = np.std(y,ddof=1)
    s = np.sqrt((n1*sig1**2 + n2*sig2**2)/(n1+n2))
    return abs(mu1 - mu2) / s

def run(_cfg,fout=None,source_data=None):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)
 
    A,C,E,W,Z,Wg = {},{},{},{},{},{}
    for i in range(1,5):
        A[i] = nx.read_graphml(cfg['refgraphs']['adj']%i)
        C[i] = nx.read_graphml(cfg['refgraphs']['chem']%i)
        E[i] = nx.read_graphml(cfg['refgraphs']['gap']%i)
        W[i] = nx.read_graphml(cfg['zhen']['white_chem']%i)
        Wg[i] = nx.read_graphml(cfg['zhen']['white_gap']%i)
        Z[i] = nx.read_graphml(cfg['zhen']['zhen_chem']%i)
    
    D = nx.DiGraph()
    for (a,b) in C[4].edges():
        if A[4].has_edge(a,b): D.add_edge(a,b)
    print(D.number_of_edges())
    print(C[4].number_of_edges(),C[4].number_of_nodes())
    #print(sorted(C[4].nodes()))
    print(A[4].number_of_edges()) 
    Ref = make_reference_graphs(C)#,remove=axon_nogo)
    print(Ref.number_of_edges())
    Ref = clean_graph(Ref,A[4])
    print(Ref.number_of_edges())
    chem = build_data(Ref,'Cook2019')
    #for c in chem: if c[1] < 0: print(c)
    #cf = pd.DataFrame(data=chem,columns=['delta','# EM sections'])

    wRef = make_reference_graphs(W)
    ref = nx.DiGraph()
    for (a,b) in wRef.edges():
        if not Ref.has_edge(a,b): continue
        ref.add_edge(a,b,weight=Ref[a][b]['weight'],id=wRef[a][b]['id'])
    
    ##Count number 1 EM sections synapses in Cook not in white
    bool_count = 0
    not_white = 0
    for (a,b) in Ref.edges():
        notwhite = not ref.has_edge(a,b)
        em1 = (Ref[a][b]['weight'] == 1)
        if notwhite: not_white += 1 
        if notwhite and em1: bool_count += 1
    print(f'Number of cook~white 1em: {bool_count}/{not_white}')

    wchem = build_data(ref,'White1986')
    #wf = pd.DataFrame(data=wchem,columns=['delta','# EM sections'])
    
    zRef = make_reference_graphs(Z)
    ref = nx.DiGraph()
    for (a,b) in wRef.edges():
        if not Ref.has_edge(a,b): continue
        ref.add_edge(a,b,weight=Ref[a][b]['weight'],id=wRef[a][b]['id'])
        
    zchem = build_data(ref,'Witvliet2020')
    zf = pd.DataFrame(data=zchem,columns=['delta','# EM sections','Data source'])


    chem = chem + wchem
    cf = pd.DataFrame(data=chem,columns=['delta','# EM sections','Data source'])
    
    wRef = make_reference_graphs(Wg)
    ref = nx.Graph()
    for (a,b) in wRef.edges():
        if not Ref.has_edge(a,b): continue
        ref.add_edge(a,b,weight=Ref[a][b]['weight'],id=wRef[a][b]['id'])
    wgap = build_data(ref,'White1986')
 
    Ref = make_reference_graphs(E)#,remove=axon_nogo)
    Ref = clean_graph(Ref,A[4])
    gap = build_data(Ref,'Cook2019')
    gap = gap + wgap
    gf = pd.DataFrame(data=gap,columns=['delta','# EM sections','Data source'])
    #gd = [cohend(gdelta[i],gdelta[j]) for (i,j) in comb]
    print(gf.tail())
    xorder = np.arange(1,5)
    corder = ['$\mathbb{C}^1$','$\mathbb{C}^2$','$\mathbb{C}^3$','$\mathbb{C}^4$']
    gorder = ['$\mathbb{G}^1$','$\mathbb{G}^2$','$\mathbb{G}^3$','$\mathbb{G}^4$']
    
    #Figure 1 
    fig,ax = plt.subplots(1,2,figsize=(14,7))
    sns.violinplot(ax=ax[0],x='delta',y='# EM sections',hue='Data source',data=cf,order=xorder,split=True)
    sns.violinplot(ax=ax[1],x='delta',y='# EM sections',hue='Data source',data=gf,order=xorder,split=True)
    ax[0].legend(loc='upper left')
    ax[0].set_xlabel(None) 
    ax[1].set_xlabel(None)
    ax[0].set_xticks(range(4))
    ax[1].set_xticks(range(4))
    ax[0].set_xticklabels(corder)
    ax[1].set_xticklabels(gorder)
    ax[0].set_ylabel('# EM sections',fontsize=18)
    ax[1].set_ylabel('# EM sections', fontsize=18)

    plt.tight_layout()
    plt.savefig('./results/syn_size_violin.png')
    
    #Figure 2 
    fig,ax  = plt.subplots(4,1,figsize=(3,4))
    cook = cf.loc[cf['Data source']=='Cook2019',['delta','# EM sections']]
    white = cf.loc[cf['Data source']=='White1986',['delta','# EM sections']]
    zhen = zf.loc[:,['delta','# EM sections']]
    for i in range(4):
        j = i + 1
        c = cook.loc[cook['delta'] == j,['# EM sections']]
        w = white.loc[white['delta'] == j,['# EM sections']]
        c.hist(ax=ax[i],bins=40,range=(0,100),grid=False,label='Cook $\mathbb{C}^%d$'%j)
        w.hist(ax=ax[i],bins=40,range=(0,100),grid=False,label='White $\mathbb{C}^%d$'%j)
        ax[i].legend(fontsize=6)
        ax[i].set_title('')
        ax[i].set_ylabel('# contacts',fontsize=7)
        ax[i].xaxis.set_minor_locator(MultipleLocator(5))
        ax[i].set_xlim([0,100])

    ax[3].set_xlabel('# EM sections on $\mathbb{M}^4$',fontsize=7)
    plt.tight_layout()
    cf.to_csv('./source_data/source_data_extended_data_f.csv',index=False)
    plt.savefig('./results/zhen_data/cook_white_histograms.svg')
    
    #Figure 3
    fig,ax = plt.subplots(1,2,figsize=(14,7))
    for i in range(4):
        j = i + 1
        c = cook.loc[cook['delta'] == j,['# EM sections']]
        z = zhen.loc[cook['delta'] == j,['# EM sections']]
        c.hist(ax=ax[0],bins=101,range=(0,101),grid=False,histtype='step',density=True,cumulative=1,
                label='$\mathbb{C}^%d$'%j,linewidth=2)
        z.hist(ax=ax[1],bins=101,range=(0,101),grid=False,histtype='step',density=True,cumulative=1,
                label='$\mathbb{C}^%d$'%j,linewidth=2)
    
    for i in range(2):
        ax[i].set_ylim([0,1])
        ax[i].set_xlim([1,100])
        ax[i].legend(loc='upper right')
        ax[i].set_xlabel('# EM sections', fontsize=12)
        ax[i].set_ylabel('Cumulative distribution',fontsize=12)
    
    ax[0].set_title('Cook2019')
    ax[1].set_title('Witvliet2020')
    plt.savefig('./results/cook_zhen_chem_em_section_survival.png')
    
    #Figure 4
    fig,ax  = plt.subplots(4,1,figsize=(3,4))
    cook = gf.loc[gf['Data source']=='Cook2019',['delta','# EM sections']]
    white = gf.loc[gf['Data source']=='White1986',['delta','# EM sections']]
    print(white.head())
    for i in range(4):
        j = i + 1
        c = cook.loc[cook['delta'] == j,['# EM sections']]
        w = white.loc[white['delta'] == j,['# EM sections']]
        c.hist(ax=ax[i],bins=40,range=(0,100),grid=False,label='Cook $\mathbb{G}^%d$'%j)
        w.hist(ax=ax[i],bins=40,range=(0,100),grid=False,label='White $\mathbb{G}^%d$'%j)
        ax[i].legend(fontsize=6)
        ax[i].set_title('')
        ax[i].set_ylabel('# contacts',fontsize=7)
        ax[i].xaxis.set_minor_locator(MultipleLocator(5))
        ax[i].set_xlim([0,100])

    ax[3].set_xlabel('# EM sections on $\mathbb{M}^4$',fontsize=7)
    plt.tight_layout()
    gf.to_csv('./source_data/source_data_extended_data_g.csv',index=False)
    plt.savefig('./results/zhen_data/cook_white_gap_histograms.svg')
    
    plt.show()

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

    run(params.config)   
