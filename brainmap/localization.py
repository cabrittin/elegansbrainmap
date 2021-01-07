"""
@name: localization.py
@description:
Module for localization functions

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from collections import defaultdict
from mpl_toolkits.axes_grid1 import make_axes_locatable


def format_loc(data,A4,nbins=51):
    m = A4.number_of_edges()
    L = np.zeros((m*nbins,4))
    bins = np.linspace(0,1,nbins)
    for (a,b,jdx,layer,sect,w) in data:
        sect = float(sect)
        if sect > 1.1 or sect < 0: continue
        jdx = int(jdx)
        idx = A4[a][b]['idx'] * nbins
        for (kdx,b) in enumerate(bins):
            if sect < b: break
        idx += kdx
        L[idx,jdx] = 1
    return L

def build_loc_data(cell,data,A4,nbins=51,return_neigh=False):
    m = A4.number_of_edges()
    L = np.zeros((m*nbins,4))
    bins = np.linspace(0,1,nbins)
    for (a,b,jdx,layer,sect,w) in data:
        sect = float(sect)
        if sect > 1.1 or sect < 0: continue
        jdx = int(jdx)
        idx = A4[a][b]['idx'] * nbins
        for (kdx,b) in enumerate(bins):
            if sect < b: break
        idx += kdx
        L[idx,jdx] = 1

    neigh = sorted([v for v in A4.neighbors(cell)])
    m = len(neigh)
    D = [np.zeros((m,nbins)) for i in range(4)]
    for (k,v) in enumerate(neigh):
        i = A4[cell][v]['idx']
        I = i * nbins
        for j in range(4):
            D[j][k,:] = L[I:I+nbins,j]
    
    if return_neigh: 
        return D,neigh 
    else:
        return D

def build_syn_data(cell,data,A4,nbins=51):
    m = A4.number_of_edges()
    L = np.zeros((m*nbins,4))
    bins = np.linspace(0,1,nbins)
    for (a,b,jdx,layer,sect,w,r) in data:
        sect = float(sect)
        if a == cell: continue
        if sect > 1.1 or sect < 0: continue
        jdx = int(jdx)
        idx = A4[a][b]['idx'] * nbins
        for (kdx,b) in enumerate(bins):
            if sect < b: break
        idx += kdx
        L[idx,jdx] = 1

    neigh = sorted([v for v in A4.neighbors(cell)])
    m = len(neigh)
    D = [np.zeros((m,nbins)) for i in range(4)]
    for (k,v) in enumerate(neigh):
        i = A4[cell][v]['idx']
        I = i * nbins
        for j in range(4):
            D[j][k,:] = L[I:I+nbins,j]
    return D
 
def raster_freq(R):
    m,n = R.shape
    rsum = R.sum(axis=0)
    f = np.zeros(m)
    for i in range(m):
        j = i + 1
        f[i] = rsum[rsum == j].sum()
    #return(f/f[0])
    return f

def raster_overlap(R):
    if R.max() == 0: return None
    r1 = max(R[0,:]*R[1,:])
    r2 = max(R[2,:]*R[3,:])
    r3 = max(R[0,:]*R[1,:]*R[2,:]*R[3,:])
    return [r1,r2,r3]

def prob_response(D):
    m,n = D[0].shape
    Z = np.zeros((m,n))
    for d in D: Z += d
    Z = Z.flatten()
    f = np.zeros(4)
    for (i,j) in enumerate(range(1,5)):
        f[i] = len(np.where(Z == j)[0])
    return f

def count_response(D):
    m,n = D[0].shape
    Z = np.zeros((m,n))
    for d in D: Z += d
    zmax = Z.max(axis=1)
    f = np.zeros(4)
    for (i,j) in enumerate(range(1,5)):
       f[i] = len(np.where(zmax == j)[0])
    return f



def row_response(D):
    m,n = D[0].shape
    k = len(D)
    F = np.zeros((m,4))
    for i in range(m):
        r = np.zeros((k,n))
        for (j,d) in enumerate(D):
            r[j,:] = d[i,:]
        F[i,:] = raster_freq(r) 
    return F

def breakdown_overlap(D):
    m,n = D[0].shape
    k = len(D)
    O = []
    for i in range(m):
        r = np.zeros((k,n))
        for (j,d) in enumerate(D):
            r[j,:] = d[i,:]
        o = raster_overlap(r)
        if not o: continue
        O.append(o)
    O = np.array(O)
    ofreq = O.sum(axis=0) / len(O)
    return ofreq

def build_lcs(X,Y):
    m = len(X)
    n = len(Y)

    L = [[None]*(n+1) for i in range(m+1)]
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0:
                L[i][j] = 0
            elif X[i-1] == Y[j-1]:
                L[i][j] = L[i-1][j-1] + 1
            else:
                L[i][j] = max(L[i-1][j], L[i][j-1])
    return L

def extract_lcs(X,Y,L):
    m = len(X)
    n = len(Y)

    index = L[m][n]
    lcs = [None]*(index + 1)
    lcs[index] = None
    
    i = m
    j = n
    while i > 0 and j > 0:
        if X[i-1] == Y[j-1]:
            lcs[index-1] = X[i-1]
            i -= 1
            j -= 1
            index -= 1
        elif L[i-1][j] > L[i][j-1]:
            i -= 1
        else:
            j -= 1
    return lcs

def format_sequences(D,neigh):
    m,n = D[0].shape
    seq = []
    for d in D:
        idx = np.where(d.transpose() > 0)
        zdx = zip(idx[0],idx[1])
        tmp = defaultdict(list)
        for (a,b) in zdx: tmp[a].append(b)
        _seq = []
        for i in sorted(tmp):
            #_seq.append(-1*i)
            _seq.extend(tmp[i])
        seq.append(_seq)
    return seq

def plot_raster(ax,D,label=None):
    for i in range(4):
        ax[i].imshow(D[i],cmap='Greys',interpolation='nearest')
        if label: ax[i].set_title(label[i],fontsize=10)
        ax[i].set_ylabel('$\mathbb{M}^4$ neighbors',fontsize=10)

def plot_contact_rate(ax,f):
    ax.bar([1,2,3,4],f,color='k',linewidth=2)
    #ax.set_ylim([0,1])
    ax.set_ylabel('Pr(# contacts $\geq x$ | 1 contact)')
    ax.set_xticks([1,2,3,4])
    ax.set_xlabel('contact rate, $x$')

def plot_contact_rate2(ax,f):
    col = ['r','#ffb311','b', 'g']
    ax.bar([1,2,3,4],f,color=col,width=0.35)
    #ax.set_ylim([0,1])
    ax.set_ylabel('# of localized contacts',fontsize=5)
    ax.set_xticks([1,2,3,4])
    ax.set_xlabel('# of datasets',fontsize=5)
 
def plot_contact_dist(ax,F):
    for i in range(1,4):
        ax.hist(F[:,i],bins=20,range=(0,1),linewidth=2,density=True,cumulative=-1,
                histtype='step',label='$\geq %d$ contacts'%(i+1))
    ax.legend(loc='lower left')
    ax.set_ylim([0,1])
    ax.set_xlim([0.01,1])
    ax.set_ylabel('Survival distribution')
    ax.set_xlabel('Contact rate broken down by neighbors')
    
def plot_overlap(ax,O): 
    ax.bar([1,2,3],O,width=0.3)
    ax.set_ylim([0,1])
    ax.set_ylabel('Fraction neighbors \nwith overlap')
    ax.set_xticks([1,2,3])
    ax.set_xlim([0.5,3.5])
    ax.set_xticklabels(['L4 L/R','Adult L/R','All'])

def plot_order(ax,lcs,size):
    Z = np.zeros(size)
    j = 0
    for l in lcs[:-1]:
        if l < 0:
            j = -1*l
            continue
        Z[l,j] = 1

    ax.imshow(Z,cmap='Greys',interpolation='nearest')
    ax.set_xlabel('Effective $z$',fontsize=7)
    ax.set_ylabel('$\mathbb{M}^4$ neighbors',fontsize=7)
    ax.set_title('LCS order',fontsize=7)


def plot_order2(ax,D,mask=None):
    #print(len(D))
    m,n = D[0].shape
    Z = np.zeros((m,n))
    for d in D: 
        Z += d
        #print(Z.min(),Z.max())
    try:
        mask[mask>0] = 1
        Z = np.multiply(Z,mask)
    except:
        pass
    
    cmap = colors.ListedColormap(['#ffffff','red', '#ffb311','blue', 'green'])
    boundaries = [-0.5,0.5,1.5,2.5,3.5,4.5]
    norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
    #print(Z.min(),Z.max())
    img = ax.imshow(Z,cmap=cmap,interpolation='nearest',norm=norm)
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    ax.set_xticks([int(j*n) for j in np.arange(0,1.2,0.2)]) 
    ax.set_xticklabels([0,0.2,0.4,0.6,0.8,1.0],fontsize=6)
    for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(6)
    cbar = plt.colorbar(img,ax=ax,fraction=0.05,pad=0.04)
    cbar.set_ticks([0,1,2,3,4]) 
    cbar.set_ticklabels(['0','1','2','3','4']) 
    cbar.ax.tick_params(labelsize=7)
    ax.set_xlabel('$\hat{z}$',fontsize=7)
    ax.set_ylabel('$\mathbb{M}^4$ neighbors',fontsize=7)
    ax.set_title('Localization of process placement',fontsize=7)
    return Z

def format_syn_data(data,Ref,rfilter=4):#thresh=4,high_pass=True):
    syn = []
    for d in data:
        if not Ref.has_edge(d[0],d[1]): continue
        #if high_pass and Ref[d[0]][d[1]]['id'] < thresh: continue
        #if not high_pass and Ref[d[0]][d[1]]['id'] > thresh: continue
        if Ref[d[0]][d[1]]['id'] != rfilter: continue 
        r = Ref[d[0]][d[1]]['id']
        syn.append(d + [r])
    return syn


