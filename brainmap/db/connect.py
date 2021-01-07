"""
db.connect.py

Module for connecting to the database.

Functions:
----------
default(db,default_file='~/.my.cnf'):
  Returns connection to database db using default setting. 
  Default settings are stored in the default file. 

connect(host,user,passwd,db):
  Returns connection to database db using passed settings.

"""
import MySQLdb


def default(db,default_file='~/.my.cnf'):
    """
    Returns connection to database db using default setting. 
    Default settings are stored in the default file.
    
    Parameters:
    ----------
    db : str
     database name
    default_file : str (default = ~/.my.cnf)
     Files with default connection settings.
    
    Returns:
    --------
    MySQLdb.connect
 
    """
    return  MySQLdb.connect(read_default_file=default_file,     
                            db=db)


def connect(host,user,passwd,db):
    """
    Returns connection to database db using passed settings.
    
    Parameters:
    -----------
    host : str
      host address
    user : str
      username
    passwd : str
      password
    db : str
      database name

    Returns:
    --------
    MySQLdb.connect    

    """

    return MySQLdb.connect(host=host, 
                           port=3306, 
                           user=user, 
                           passwd=passwd, 
                           db=db)
    
