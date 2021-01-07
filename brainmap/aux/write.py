"""
write.py

Part of auxilary module. Contains functions for writing data.

Required 3rd party packages:
  csv

Author: Christopher Brittin

"""

import csv

def from_dict(fout,Dic):
    """
    Write data from dictionary to file
    
    Parameters
    -----------
    fout : str
       Path to output file
    Dic : dictionary
       Dictionary to be written to data. Treated as 1D dictionary, meaning 
       each key has one value. 

    """
    fout = open(fout,'w')
    fWriter = csv.writer(fout,delimiter = ',',
                         quotechar = '"',
                         quoting = csv.QUOTE_MINIMAL)
    for k in sorted(Dic.keys()):
        if type(Dic[k]) == list:
            fWriter.writerow([k] + Dic[k])
        else:
            fWriter.writerow([k,Dic[k]])

def from_list(fOut,List):
    """
    Write data from list to file
    
    Parameters
    -----------
    fOut : str
       Path to output file
    List : list
       List to be written to data. Each row can have multiple elements. 

    """

    fOut = open(fOut,'w')
    if type(List[0]) in (list,tuple):
        fWriter = csv.writer(fOut,delimiter = ',',
                             quotechar = '"',
                             quoting = csv.QUOTE_MINIMAL)
        for l in List:
            fWriter.writerow(l)
    else:
        for l in List:
            fOut.write(''.join([str(l),'\n']))
     
    fOut.close()
