"""
@name: segment_analysis.py
@description:
Module for analysis of bundle contact at EM segments

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05
"""

from lxml import etree
import numpy as np
from tqdm import tqdm
from collections import defaultdict

import aux

def extract_bundle_contact(G,xml,rlmap):
    num_class = len(set([G.node[n]['bid'] for n in G.nodes()]))
    tree = etree.parse(xml)
    root = tree.getroot()
    cells = sorted([c.get('name') for c in root.findall('cell')])
    data = defaultdict(list)
    for c in tqdm(cells,desc='Cells'):
        _c = c
        if _c in rlmap: _c = rlmap[_c]
        if not G.has_node(_c):continue
        _data = extract_cell_bundle_contact(c,G,root,num_class,rlmap)
        data[_c] += _data

    return data

def extract_cell_bundle_contact(c,G,root,num_class,rlmap):
    _c = c
    if _c in rlmap: _c = rlmap[_c]
    cls = G.node[_c]['bid']
    cell = root.find("cell[@name='%s']" %c)
    data = []
    for layer in cell.findall('layer'):
        for idx in layer.findall('idx'):
            z = np.zeros(num_class)
            num_neigh = 0
            for n in idx.findall('neighbor'):
                neigh = n.get('name')
                _neigh = neigh
                if _neigh in rlmap: _neigh = rlmap[_neigh]
                if not G.has_edge(_c,_neigh): continue
                _cls = G.node[_neigh]['bid']
                adj = int(n.get('adjacency'))
                z[_cls] += adj 
                num_neigh += 1
                #print(c,layer.get('name'),idx.get('value'),neigh,adj)
            #if cls < 5: z[5] = 0
            if z.sum() < 1: continue 
            z /= z.sum()
            #tmp = [c,cls,z[cls],num_neigh,layer.get('name'),idx.get('value')] + z.tolist()
            tmp = [c,layer.get('name'),idx.get('value'),num_neigh] + z.tolist()
            data.append(tmp)

    return data


def get_cell_contact(flist):
    z =[]
    for f in flist: z += aux.read.into_list2(f)
    n = len(z[0]) - 4
    m = len(z)
    Z = np.zeros((m,n))
    for (i,l) in enumerate(z): Z[i,:] = l[4:] 
    return Z 

