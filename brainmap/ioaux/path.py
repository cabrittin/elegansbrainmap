"""
path

submodule for traversing file system
"""

import os

def ls(path):
    """
    Iterate through files in a directory
    """
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path,file)):
            yield file
