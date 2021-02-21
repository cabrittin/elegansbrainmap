"""
tpair_adj_deg.py

Plots distribution of the differences in homologous adjacency degrees.

created: Christopher Brittin
date: 01 November 2018

"""
import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind
from scipy.stats import wilcoxon
from scipy.stats import ttest_1samp
import pandas as pd
import seaborn as sns
import matplotlib as mpl

#Brittin modules
from connectome.load import from_db
from networks.stats import *
from plots import print_wilcoxon,boxplots
#from figures.stats import *
import aux

mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5


CONFIG = os.environ['CONFIG']
REMOVE = ['CEHDL','CEHDR','CEHVL','CEHVR','HSNL','HSNR','PLNL','PLNR','PVNL','PVNR','PVPL','PVPR']
ADULT_COL = '#FFC300'
L4_COL = '#3380FF'
AL_COL = '#14FA29'

def write_out(nodes,score1,score2,_fout):
    _out = []
    for (n1,n2) in nodes:
        _out.append([n1,n2,score1[n1],score2[n2],abs(score1[n1]-score2[n2])])
    aux.write.from_list(_fout,_out)

def run(_cfg,fout=None,source_data=None):
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(_cfg)
    
    #_remove = aux.read.into_list(cfg['mat']['remove']) 
    _remove = ['VC01','VD01','VB01','VB02','HSNL','HSNR','PVNL','PVNR','PLNL','PLNR','PVR','PVR.']
    left =  aux.read.into_list(cfg['mat']['left_nodes']) 
    right = aux.read.into_list(cfg['mat']['right_nodes'])
    lrmap = aux.read.into_lr_dict(cfg['mat']['lrmap']) 
    data = []

    N2U = 'N2U'
    JSH = 'JSH'
    n2u = from_db(N2U,adjacency=True,remove=_remove)
    jsh = from_db(JSH,adjacency=True,remove=_remove)
    ndelta,jdelta,bdelta = [],[],[]

    lnd = get_adj_deg(n2u,vertices = left)
    rnd = get_adj_deg(n2u,vertices = right)
    tmp = [n for n in sorted(lnd)]
    for (l,r) in [(n,lrmap[n]) for n in sorted(lnd.keys())]: 
        data.append(['Adult L/R',l,r,lnd[l],rnd[r],lnd[l]-rnd[r]])
        ndelta.append(lnd[l]-rnd[r])
    
    lnd = get_adj_deg(jsh,vertices = left)
    rnd = get_adj_deg(jsh,vertices = right)
    for (l,r) in [(n,lrmap[n]) for n in sorted(lnd.keys())]: 
        data.append(['L4 L/R',l,r,lnd[l],rnd[r],lnd[l]-rnd[r]])
        jdelta.append(lnd[l]-rnd[r])
    
    cells = []
    for n in sorted(lnd.keys()):
        cells.append(n)
        cells.append(lrmap[n])
    bnd = get_adj_deg(n2u,vertices = cells)
    bjd = get_adj_deg(jsh,vertices = cells)
    for c in cells:
        data.append(['Adult/L4',c,c,bnd[c],bjd[c],bnd[c]-bjd[c]])
        bdelta.append(bnd[c]-bjd[c])
     
    df = pd.DataFrame(data,columns=["Comparison","Cell1","Cell2","Deg1","Deg2","Deg_diff"])
    print('Stats:')
    print_wilcoxon(ndelta,'Adult L/R')
    print_wilcoxon(jdelta,'L4 L/R')
    print_wilcoxon(bdelta,'Adult/L4',alternative="greater")
    
    #tval1,pval1 = ttest_ind(ndelta,jdelta)
    #tval2,pval2 = ttest_ind(jdelta,bdelta)
    #tval3,pval3 = ttest_ind(ndelta,bdelta)
    
    sns.set_theme(style="whitegrid")
    fig,ax = plt.subplots(1,1,figsize=(2.15,1.7))
    flierprops = dict(markersize=1,marker='d',markerfacecolor='k')
    medianprops = dict(linestyle='-',linewidth=0.5,color='k')
    whiskerprops = dict(linestyle='-',linewidth=0.3,color='k')
    capprops = dict(linewidth=0.3)
    sns.boxplot(x="Comparison",y="Deg_diff",
            data=df,width=0.3,ax=ax,linewidth=0.3,color="#a5a5a5",
            flierprops=flierprops,medianprops=medianprops,capprops=capprops)
    ax.set_ylim([-30,30])
    ax.set_yticks([-30,-20,-10,0,10,20,30])
    #ax.set_yticklabels([-30,-20,-10,0,10,20,30],fontsize=5)
    for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(7)
    for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(5)
    ax.axhline(0,color='r',linewidth=0.8,linestyle='--')
    ax.set_xlabel("")
    ax.set_ylabel("Degree difference",fontsize=7)
    plt.tight_layout() 
    if fout: plt.savefig(fout)
    plt.show()
    if source_data: df.to_csv(source_data,index=False)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    
    parser.add_argument('-o','--fout',
                        action = 'store',
                        dest = 'fout',
                        default = None,
                        required = False,
                        help = 'Output svg file')
    
    parser.add_argument('-sd','--source_data',
                        action = 'store',
                        dest = 'source_data',
                        default = None,
                        required = False,
                        help = 'Output source data file')
    
    params = parser.parse_args()
    run(params.config,fout=params.fout,source_data = params.source_data)
