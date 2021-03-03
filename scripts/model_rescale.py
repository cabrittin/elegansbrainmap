"""
@name: prec_rescale.py
@description:
    Rescale A1, A2, A3


@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-03
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

FIN = './data/model/model_%d.npy'
FOUT = './results/rescaled_contacts.svg'
P = np.array([0.93,0.91])
S = np.array([0.72,0.81])
F = np.array([0.33,0.23])
#VAR = 4*(F*P*(1-P) + (1-F)*S*(1-S))
VAR = np.array([1,1])

mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 5

def build_data(M,idx):
    D = np.zeros((4,4,))
    for (i,j) in enumerate(range(4,0,-1)):
        D[i,:] = M[j][idx]
    D[1,3:] = 0
    D[2,2:] = 0
    D[3,1:] = 0
    return D

def predicted(D):
    scale = D.sum(axis=1)
    print(D)
    print('scale',scale)
    print('scale',scale[:,None],scale[0])
    R = np.tile(D[0,:],(4,1)) * scale[:,None] / scale[0]
    R[D==0] = 0
    return R 

def plot_rescale(ax,A,R,err,width=0.25,width2=0.23,labels=None,ylabel=[None,None,None]):
    aerr = np.sqrt(A*err)
    rerr = np.sqrt(R*err)

    for (j,i) in enumerate(range(1,4)):
        #print(f'axis: {j}')
        #print('actual',A[i,:])
        #print('predicted',R[i,:])
        #print('actual_err',aerr)
        #print('pred_err',rerr)
        _plot_rescale(ax[j],A[i,:],R[i,:],aerr[i,:],rerr[i,:],labels=labels,ylabel=ylabel[j])

def _plot_rescale(ax,data,pred,err1,err2,width=0.25,width2=0.23,labels=None,ylabel='#contacts'):
    ind = np.arange(4)
    #ax.grid(zorder=-1)
    ax.bar(ind-0.5*width,data,yerr=err1,width=width,ecolor='black', capsize=1,label='Empirical')
    ax.bar(ind+0.5*width,pred,yerr=err2,width=width2,ecolor='black', capsize=1,label='Predicted')
    
    ax.set_xlim([-0.5,3.5])
    ax.set_xticks(ind)
    if labels: ax.set_xticklabels(labels)
    ax.legend(fontsize=6)
    ax.set_ylabel('# contacts',fontsize=7) 
    ax.set_ylabel('# contacts',fontsize=7) 
    if ylabel: ax.set_title(ylabel,fontsize=7)
    
def run(cfg,fout=None,source_data=None):
    
    M = dict([(i,np.load(FIN%i)) for i in range(1,5)])

    print(VAR)
    
    #print('M3',M[3])
    #print('M4',M[4])
    for i in range(1,5): print(f'M{i}',M[i])
    C = build_data(M,1)
    RC = predicted(C) 
    print('RC',RC)

    G = build_data(M,2)
    RG = predicted(G)

    _label = ['$\mathbb{%s}^1$','$\mathbb{%s}^2$','$\mathbb{%s}^3$','$\mathbb{%s}^4$']
    clabel = [l%'C' for l in _label]
    glabel = [l%'G' for l in _label]

    #s = '# contacts scaled by %s' 
    s = 'scaled by %s' 
    ylabel = [s%'$\mathbb{M}^3$',s%'$\mathbb{M}^2$',s%'$\mathbb{M}^1$']

    fig,ax = plt.subplots(3,2,figsize=(3.5,4))
    plot_rescale(ax[:,0],C,RC,VAR[0],labels=clabel,ylabel=ylabel)
    plot_rescale(ax[:,1],G,RG,VAR[1],labels=glabel,ylabel=ylabel)
    #ax[0,0].set_title('Chemical synapses',fontsize=7)
    #ax[0,1].set_title('Gap junctions',fontsize=7)
    plt.tight_layout()
    plt.savefig(FOUT) 
    plt.show()

if __name__=="__main__":
    run(1)

