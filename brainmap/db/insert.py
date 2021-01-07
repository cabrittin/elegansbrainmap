"""
insert.py

Module for inserting data into database

Functions:
---------
adjacency(cur,data)
   Inserts adjacency data

radial_pharynx(cur,data):
   Insert radial pharynx data [OBJ_Name,radial_distance,phi]

created: Christopher Brittin
date: 01 November 2018
"""

def adjacency(cur,data):
    """
    Inserts adjacency data

    Parameters:
    -----------
    cur : MySQLdb cursors
    data : list
      Adjacency data. Row format
      [cell1,cell2,idx1,idx2,imgNum,adjacency,obj1,obj2]
    
    """
    sql = ("insert into adjacency2 "
           "(pre,post,preidx,postidx,imgNum,weight,preObj,postObj) "
           "values (%s,%s,%s,%s,%s,%s,%s,%s)")
    cur.executemany(sql,data)
    
def radial_pharynx(cur,data):
    """
    Insert radial pharynx data

    Parameters:
    -----------
    cur : MySQLdb cursors
    data : list
      Adjacency data. Row format
       [OBJ_Name,radial_distance,phi]    
    """
    sql = ("insert into radialPharynx "
           "(Obj_Name,distance,phi) "
           "values (%s,%s,%s)")    
    cur.executemany(sql,data)
