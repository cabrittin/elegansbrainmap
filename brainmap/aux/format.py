"""
dtypes.py
Author: Christopher Brirint

Functions for manipulating datatypes

"""


import re

def rm_brack(s):
    #Removes surrounding brackets of a
    #string e.g. '[a]' -> 'a' 
    return re.sub('[\[\]]','',s)

def get_args(kwargs,default):
    #Sorts keyword args into dict
    #kwargs: keyword arguments
    #default: dictionary of default values
    args = {}
    for d in default:
        if d in kwargs and kwargs[d]:
            args[d] = (kwargs[d])
        else:
            args[d] = (default[d])
    return args

def chunk_list(l,n):
    #Chunks list l into smaller lists
    #of size n
    for i in range(0,len(l),n):
        yield l[i:i+n]

def swap_var(x,y):
    #swaps x for y
    temp = x
    x = y
    y = temp
    return x,y

def sort_dictionary_by_value(D,reverse):
    #Sorts dictionary keys by value
    import operator
    return sorted(D.iteritems(), key=operator.itemgetter(1),reverse=reverse)

def merge_two_dicts(a, b):
    '''Given two dicts, merge them into a new dict'''
    r = dict(a.items() + b.items() +
             [(k, a[k] + b[k]) for k in set(b) & set(a)])
    return r 

def switch_keys_vals(ncat):
    nlist = {}
    for n in ncat:
        c = ncat[n]
        if c not in nlist: nlist[c] = []
        nlist[c].append(n)
    return nlist

class ScaleRange:
    def __init__(self,xmin,xmax,ymin,ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.m = float(ymax-ymin)/(xmax-xmin)
        
    def linear_scale(self,x):
        f = self.m*(x-self.xmin) + self.ymin
        return f
