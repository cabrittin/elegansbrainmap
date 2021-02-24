"""
@name: expand_list_to_bilateral.py
@description:
Expand list with only left cells to include both left and right cells

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-05
"""

import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation

import aux


CONFIG = os.environ['CONFIG']

def expand_list(lst,lrdict):
    new = []
    for [a,b] in lst:
        new.append([a,b])
        ar = lrdict[a]
        if a == ar: continue
        new.append([ar,b])
    return new



if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('lst',
                        action = 'store',
                        help = 'List to expand')

    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    
    params = parser.parse_args()
    
    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)

    lrdict = aux.read.into_lr_dict(cfg['mat']['lrmap'])
    
    fout = params.lst.replace('.','_lr.')
    lst = aux.read.into_list2(params.lst)
    new = expand_list(lst,lrdict)
    aux.write.from_list(fout,new)
