"""
paper_figure.py

Generates plots and tables for:

Brittin, C. A. , Cook, S. J., Hall, D.H., Emmons, S. W., Cohen. N. A multiscale brain map derived from whole-brain volumetric reconstructions. Nature (2021).

created: Christopher Brittin
date: February 2021 

Synopsis:
  python paper_figures figNum [fout]


Parameters:
  fig (str): Figure number or table number from manuscript.
             e.g. f1a is Figure 1a
                  fs2a is Figure S2a
                  ts1 is Table S1
  -o, --output (str): Output file name (optional).
  -s, --source_data (str): Output file name for source data (optional)

"""

import os
import sys
sys.path.append(r'./scripts')
import argparse
import importlib

import aux

#CONFIG = os.environ['CONFIG']
CONFIG = 'configs/config.ini'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fig',
                        action="store",
                        help=("Figure number or table number."
                              "Examples: "
                              "To generate Figure 1a pass f1a. "
                              "To generate Extended Data 2a pass ed2a. "
                              "To generate Table S1 pass ts1." )
                        )
    
    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
   
    parser.add_argument('-o','--output',
                        dest = 'fout',
                        action="store",
                        required= False,
                        default = None,
                        help="Output graph file"
                        )
    
    parser.add_argument('-s','--source_data',
                        dest = 'source_data',
                        action="store",
                        required= False,
                        default = None,
                        help="Output source data file."
                        )
    

    params = parser.parse_args()

    figmap = aux.read.into_dict('paper_figmap.csv')
    
    print(params.config)

    if params.fig in figmap:
        module = importlib.import_module(figmap[params.fig]) 
        module.run(params.config,fout=params.fout,source_data = params.source_data)

    else:
        print("Error: Not a valid paper figure.")
