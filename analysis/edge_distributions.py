"""
edge_distributions.py

Generate plots concerning the edge distributions

crated: Christopher Brittin
data: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric')
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes, InsetPosition
import matplotlib as mpl
import numpy as np

#Brittin modules
from connectome.load import from_db
from plots import plot_dist
import ioaux

mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5


REMOVE = ['VC01','VD01','VB01','VB02']
SCALE = 5*90*(1e-6)
LOWCOL = '#42fc30'
HIGHCOL = '#7c3aff'
MIDCOL = '#bebebe'

def run(db='N2U',fout='results/N2U_edge_dist.svg',source_data=None):
    C = from_db(db,adjacency=True,remove=REMOVE)
    edges = [e['weight']*SCALE for e in C.A.es()]
    _edges = [[C.A.vs[e.source]['name'],C.A.vs[e.target]['name'],e['weight']*SCALE] for e in C.A.es()]
    _edges = sorted(_edges)
    print(_edges[0],_edges[-1])
    _edges = [['cell1','cell2','membrane_contact_microns^2']] + _edges

    print(np.max(edges))
    low = np.percentile(edges,35)
    high = np.percentile(edges,66)
    print(low,high)
    print('# edges',len(edges))

    fig,ax = plt.subplots(1,1,figsize=(2,1.8))
    plot_dist(ax,edges,density=True,cumulative=-1,xlim=[0,10],ylim=[0,1],
              hrange=(0,60),nbins=1000,
              xlabel='Membrane contact area ($\mu$m$^2$)',
              ylabel='Survival distribution',fs=7)
    hist, bin_edges = np.histogram(edges,bins=1000,range=(0,60),density=True,normed=True)
    dx = bin_edges[1] - bin_edges[0]
    hist = 1 - np.cumsum(hist*dx)
    ax.fill_between(bin_edges[:-1],hist,color=MIDCOL)
    ax.fill_between(bin_edges[:-1],hist, where = hist < 0.35,color=HIGHCOL)
    ax.fill_between(bin_edges[:-1],hist, where = hist > 0.66,color=LOWCOL)
    

    edges = np.array(edges)
    esum = np.sum(edges)
    elow = np.sum(edges[np.where(edges < low)]) / esum
    ehigh = np.sum(edges[np.where(edges > high)]) / esum
    emid = 1 - elow - ehigh 
    #fig,ax = plt.subplots(1,1)
    ax2 = plt.axes([0,0,1,1])
    ip = InsetPosition(ax, [0.1,0.4,0.5,0.5])
    ax2.set_axes_locator(ip)
    wedges, texts, autotexts = ax2.pie([elow,emid,ehigh],textprops={'fontsize':6},
            colors=[LOWCOL,MIDCOL,HIGHCOL],autopct='%1.1f%%')
    legend = ax2.legend(wedges,['Low','Mid','High'],title='Contact range',
                loc='center left',bbox_to_anchor=(1, 0, 0.5, 1),fontsize=6)
    ax2.set_title('Total surface area',fontsize=6)
    #plt.setp(autotexts, size=12, weight="bold")
    plt.setp(legend.get_title(),fontsize=6)
    plt.tight_layout()
    if fout: plt.savefig(fout)    
    plt.show()
    if source_data: ioaux.write.from_list(source_data,_edges)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('db',
                        action="store",
                        help=("Database name")
                        )
   
    parser.add_argument('-o','--output',
                        dest = 'fout',
                        action="store",
                        required= False,
                        default = None,
                        help="Output file, if you wish to save the output."
                        )
    
    parser.add_argument('-s','--source_data',
                        dest = 'source_data',
                        action="store",
                        required= False,
                        default = None,
                        help="Source data file"
                        )
    
    params = parser.parse_args()

    run(db=params.db,fout=params.fout,source_data=params.source_data)


