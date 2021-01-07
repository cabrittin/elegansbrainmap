"""
read.py

Part of auxilary module. Contains functions for reading data from files.

Required 3rd party packages:
  csv

Author: Christopher Brittin

"""

import csv

def into_dict(fIn,delimiter=","):
    """
    Read data from file into a dictionary. By default
    first element of each row is assigned to the key
    and the second element is assigned to the value.

    Parameters
    ----------
    fIn : str
      Path to input file
    delimeter : str (optional)
      Delimeter for parsing lines in file. (default is ',')
    
    """
    fIn = open(fIn,'r')
    fInReader = csv.reader(
        fIn,delimiter = delimiter, 
        quotechar = ' ', quoting = csv.QUOTE_NONE)
    d = {}
    for line in fInReader:
        if line: d[line[0]] = line[1]
    fIn.close()
    return d  
  
def into_dict2(fIn,delimiter=","):
    """
    Read data from file into a 2D dictionary. By default
    first element of each row is assigned to the key
    and the remaining elements are made into a list
    and assigned to the value.

    Parameters
    ----------
    fIn : str
      Path to input file
    delimeter : str (optional)
      Delimeter for parsing lines in file. (default is ',')
    
    """
    
    fIn = open(fIn,'r')
    fInReader = csv.reader(
        fIn,delimiter = delimiter, 
        quotechar = ' ', quoting = csv.QUOTE_NONE)    
    d = {}
    for line in fInReader:
        if line: d[line[0]] = line[1:]
    fIn.close()
    return d 


def into_list(fin):
    """
    Read data from file into a list.
    Lines with '#' are not read.

    Parameters
    ----------
    fIn : str
      Path to input file
     """
    fIn = open(fin,'r')
    fInReader = csv.reader(
        fIn,delimiter = ',', quotechar = ' ', quoting = csv.QUOTE_NONE)
    lst = [l[0] for l in fInReader if '#' not in l]
    fIn.close()
    return lst

def into_list2(fIn,**kwargs):
    """
    Read data from file into a 2D list.
    Lines with '#' are not read.

    Parameters
    ----------
    fIn : str
        Path to input file
    delimeter : str (optional)
        Delimeter for parsing lines in file. (default is ',')
    type : str (optional)
        Converst list to specified data type. 
        Choices 'int', 'float' and 'str'. 
    """   
    delimiter = ','
    if 'delimiter' in kwargs:
        delimiter = kwargs['delimiter']
    tmap = {'int':int,'float':float,'str':str}
    fIn = open(fIn,'r')
    fInReader = csv.reader(
        fIn,delimiter = delimiter, quotechar = ' ', quoting = csv.QUOTE_NONE)
    if 'type' in kwargs and kwargs['type'] in kwargs.values():
        return [map(tmap[kwargs['type']],row) for row in fInReader if '#' not in row[0]] 
    else:
        return [row for row in fInReader if '#' not in row[0] ]

  
def into_lr_dict(fin):
    """
    Creates left/rigth dictionary for cells

    Parameters
    ----------
    fIn : str
      Path to file specifying left/right cells. Should have format
      'left_cell,right_cell'
    """   
    lr = into_dict(fin)
    _keys = list(lr.keys())
    for key in _keys:
        lr[lr[key]] = key
    return lr

def into_homolog_dict(fin):
    """
    Creates homologous dictionary for cells

    Parameters
    ----------
    fIn : str
      Path to file specifying homologous cells. Should have format
      'homolog_class,cell_name1,cell_name2,...'
    """       
    nc = into_list2(fin)
    nclass = {}
    for row in nc:
        if len(row) == 1:
            nclass[row[0]] = row[0]
        else:
            for key in row[1:]:
                nclass[key] = row[0]
    return nclass    


def into_map(fIn,**kwargs):
    """
    Creates a dictionary map where all elements in the row are mapped to 
    the first element in the row

    Parameters
    ----------
    fIn : str
      path to specify map file
    delimeter : str (optional)
        Delimeter for parsing lines in file. (default is ',')
    """       
    delimiter = ','
    if 'delimiter' in kwargs:
        delimiter = kwargs['delimiter']

    fIn = open(fIn,'r')
    fInReader = csv.reader(
        fIn,delimiter = delimiter, 
        quotechar = ' ', quoting = csv.QUOTE_NONE)    
    d = {}
    for line in fInReader:
        for l in line[1:]:
            d[l] = line[0]
    fIn.close()
    return d     
