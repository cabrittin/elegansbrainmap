"""
@name:
@description:
Pull A4 adjacency data

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05
"""

import os
from configparser import ConfigParser,ExtendedInterpolation
import argparse
import networkx as nx
from tqdm import tqdm


import db as db
from connectome.load import from_db
import ioaux
from connectome.format_graphs import filter_graph_edge

#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'

def format_img(img,_db):
    if _db == 'JSH':
        dimg = dict([(im[3:],i) for (i,im) in  enumerate(img)])
    elif _db == 'N2U':
        dimg = {}
        for (i,im) in enumerate(img):
            if 'NR' in im: im = im.replace('NR','_')
            if 'VC' in im: im = im.replace('VC','_VC_')
            dimg[im] = i
        dimg['N2U_0164'] = dimg['N2U_164']
    return dimg

def pull_adjacency(data,_db,A4,_min=0,_max=500,imap={'JSH':[0,1],'N2U':[2,3]}):
    D = from_db(_db,adjacency=True,remove=remove,dataType='networkx')
    D.A = filter_graph_edge(D.A,pct=edge_thresh)
    D.split_left_right(left,right)  
    #D.map_right_graphs(lrmap)
   
    con = db.connect.default(_db)
    cur = con.cursor()
    adj = db.mine.get_adjacency_data(cur)
    
    img = db.mine.get_img_number(cur)
    dimg = format_img(img,_db)    
    tst = ['AVKR','SIAVR']
    for (a,b,w,layer) in tqdm(adj,desc="Adjacency iter:"):
        sect = dimg[layer]
        if sect < _min or sect > _max: continue
        ra = lrmap[a]
        rb = lrmap[b]
        inleft = D.Al.has_edge(a,b) and A4.has_edge(a,b)
        inright = D.Ar.has_edge(a,b) and A4.has_edge(ra,rb)
        if inleft: 
            idx = imap[_db][0]
            _sect = section_map(idx,sect)
            data.append([a,b,idx,layer,_sect,w])
        if inright: 
            idx = imap[_db][1]
            _sect = section_map(idx,sect)
            data.append([ra,rb,idx,layer,_sect,w])
        #if a in tst and b in tst:
        #    print(a,b,w,layer,dimg[layer],inleft,inright,ra,rb,A4.has_edge(ra,rb),D.Ar.has_edge(a,b))
 
        #print(layer,dimg[layer])

def pull_synapses(data,_db,A4,_min=0,_max=500,imap={'JSH':[0,1],'N2U':[2,3]}):
    D = from_db(_db,chemical=True,remove=remove,dataType='networkx')
    D.split_left_right(left,right)  
    #D.map_right_graphs(lrmap)
   
    con = db.connect.default(_db)
    cur = con.cursor()
    adj = db.mine.get_synapse_data(cur,'chemical')
    
    for (a,post,sections,continNum,series) in tqdm(adj,desc="Chemical iter"):
        if a not in lrmap: continue 
        for (oname,x,y,sect,ifile) in db.mine.get_contin_xyz(cur,continNum):
            if sect < _min or sect > _max: continue
            for b in post.split(','):
                if b not in lrmap: continue
                ra = lrmap[a]
                rb = lrmap[b]
                inleft = D.Cl.has_edge(a,b) and A4.has_edge(a,b)
                inright = D.Cr.has_edge(a,b) and A4.has_edge(ra,rb)
                if inleft: 
                    idx = imap[_db][0]
                    _sect = section_map(idx,sect)
                    data.append([a,b,idx,ifile.split('.')[0],_sect,continNum])
                if inright: 
                    idx = imap[_db][1]
                    _sect = section_map(idx,sect)
                    data.append([ra,rb,idx,ifile.split('.')[0],_sect,continNum])

def pull_gap_junctions(data,_db,A4,_min=0,_max=500,imap={'JSH':[0,1],'N2U':[2,3]}):
    D = from_db(_db,electrical=True,remove=remove,dataType='networkx')
    D.split_left_right(left,right)  
    #D.map_right_graphs(lrmap)
   
    con = db.connect.default(_db)
    cur = con.cursor()
    adj = db.mine.get_synapse_data(cur,'electrical')
     
    for (a,b,sections,continNum,series) in tqdm(adj,desc='Gap j. iter:'):
        if a not in lrmap: continue 
        if b not in lrmap: continue 
        for (oname,x,y,sect,ifile) in db.mine.get_contin_xyz(cur,continNum):
            if sect < _min or sect > _max: continue
            ra = lrmap[a]
            rb = lrmap[b]
            inleft = D.El.has_edge(a,b) and A4.has_edge(a,b)
            inright = D.Er.has_edge(a,b) and A4.has_edge(ra,rb)
            if inleft: 
                idx = imap[_db][0]
                _sect = section_map(idx,sect)
                data.append([a,b,idx,ifile.split('.')[0],_sect,continNum])
            if inright: 
                idx = imap[_db][1]
                _sect = section_map(idx,sect)
                data.append([ra,rb,idx,ifile.split('.')[0],_sect,continNum])


def section_map(idx,_sect):
    if idx == 0:
        sect = map_to_interval(_sect)
    elif idx == 1:
        sect1 = map_l4_right_to_left(_sect)
        sect = map_to_interval(sect1)
    elif idx == 2:
        sect1 = map_adult_left_to_l4_left(_sect)
        sect = map_to_interval(sect1)
        #if _sect < 186:
        #    print(_sect,sect1,sect)
    elif idx == 3:
        sect1 = map_adult_right_to_l4_left(_sect)
        sect = map_to_interval(sect1)
        #if _sect < 186:
        #    print(_sect,sect1,sect)
    return sect

def map_to_interval(x):
    m = 1. / 280
    b = 15*m
    y = m * x - b
    return y

def map_l4_right_to_left(x):
    y = x
    if x < 172:
        m = (158. - 97) / (171. - 108)
        b = 97 - m * 108 
        y = m * x + b
    elif x < 190:
        m = (190.-158) / (190 - 172)
        b = 158 - m*172
        y = m * x + b
    return y

def map_adult_left_to_l4_left(x):
    y = x
    if x < 157:
        m = (158. - 97) / (147. - 94)
        b = 97 - m * 94
        y = m * x + b
    elif x < 182:
        m = (174. - 158) / (157. - 147)
        b = 158 - m * 147
        y = m * x + b
    else:
        m = (295. - 214) / (230. - 186)
        b = 214 - m * 186
        y = m * x + b
    return y

def map_adult_right_to_l4_left(x):
    y = x
    if x < 157:
        m = (158. - 97) / (133. - 87)
        b = 97 - m * 87
        y = m * x + b
    elif x < 182:
        m = (174. - 158) / (157. - 147)
        b = 158 - m * 147
        y = m * x + b
    else:
        m = (295. - 214) / (230. - 186)
        b = 214 - m * 186
        y = m * x + b
    return y


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
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
    left = ioaux.read.into_list(cfg['mat']['left_nodes'])
    right = ioaux.read.into_list(cfg['mat']['right_nodes'])
    lrmap = ioaux.read.into_lr_dict(cfg['mat']['lrmap'])
    remove = ioaux.read.into_list(cfg['mat']['remove'])
    edge_thresh = cfg.getint('params','lower_weight_threshold')
    
    A4 = nx.read_graphml(cfg['refgraphs']['adj']%4)
    
    data = []
    chem = []
    gap = []

    _db = 'JSH'
    _min = cfg.getint('adj_align','%s_min'%_db)
    _max = cfg.getint('adj_align','%s_max'%_db)
    pull_adjacency(data,_db,A4,_min=_min,_max=_max)
    pull_synapses(chem,_db,A4,_min=_min,_max=_max)
    pull_gap_junctions(gap,_db,A4,_min=_min,_max=_max) 
    _db = 'N2U'
    _min = cfg.getint('adj_align','%s_min'%_db)
    _max = cfg.getint('adj_align','%s_max'%_db)
    pull_adjacency(data,_db,A4,_min=_min,_max=_max)
    pull_synapses(chem,_db,A4,_min=_min,_max=_max)
    pull_gap_junctions(gap,_db,A4,_min=_min,_max=_max) 

    ioaux.write.from_list(cfg['adj_align']['fout'],data)
    ioaux.write.from_list(cfg['adj_align']['cout'],chem)
    ioaux.write.from_list(cfg['adj_align']['gout'],gap)
    
    

