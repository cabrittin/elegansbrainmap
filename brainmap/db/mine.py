"""
mine.py

Module for extracting and writing data to the database

Functions:
---------
get_db(cur)
   Returns database name

get_neurons(cur)
   Returns list of all cell names

get_contins(cur,cell)
   Returns list of contins for cell name

get_img_number(cur,**kwargs)
   Returns list of image(sections) names from the NR an VC
   series. Use kwargs to define range = (start,end) where 
   start and end are ints.

get_obj_area(cur)
   Returns list of objects and associated areas

get_synapse_contins(cur,stype,**kwargs)
   Returns list of synapse contin numbers

get_synapse_data(cur,stype,**kwargs)
   Returns list of synapses

get_adjacency_data(cur)
   Returns list of adjacency data

get_adjacency_cells(cur)
   Return list of cell in the adjacency data

get_touch_density(cur,key="pixels_norm")
   Returns list of 'touch_density' measures. Possible key
   values are [pixels,pixels_norm,segments,segments_norm]

maxAnterior(cur,cell,contin=None)
   Returns the max anterior section of cell

maxPosterior(cur,cell,contin=None)
   Returns the max posterior sections of cell

neuron_cylinder(cur)
   Returns the cylindrical locations of neuron segment centroids

gap_cylinder(cur)
   Return the cylindrical coordinates of the gap junction

syn_cylinder(cur)
   Return the cylindrical coordinates of synapses

get_cell_genes(cur,cell):
    Return the genes expressed in cell

get_object_xy_in_layer(cur,cell,layer):
    Returns (object_num,x,y) coordinates of cell object in layer

get_synapse_from_layer(cur,layer):
    Returns synapses in given layer

get_object_adjacency(cur,obj):
    Returns the adjacency of the object number

order_synapses_by_section_number(cur,stype='chemical')
    Returns list of synapses ordered by sections number 

created: Christopher Brittin
date: 01 November 2018
"""


import aux

def get_db(cur):
    """
    Returns database name

    Parameters:
    -----------
    cur : MySQLdb cursors

    """
    cur.execute("select database()")
    return cur.fetchone()[0]

def get_neurons(cur):
    """
    Returns list of all contin names

    Parameters:
    -----------
    cur : MySQLdb cursors

    """
    sql = ("select distinct(CON_AlternateName) "
           "from contin "
           "where type like 'neuron' "
           "and CON_Remarks like '%OK%'")
    cur.execute(sql)
    
    return list(set([aux.format.rm_brack(a[0]) for a in cur.fetchall()]))

def get_contins(cur,cell):
    """
    Returns list of all contin names

    Parameters:
    -----------
    cur : MySQLdb cursors
    cell : str
      Cell name
    """
    sql = ("select CON_Number "
           "from contin "
           "where CON_AlternateName like '%s'"
           %cell)
    cur.execute(sql)
    return [a[0] for a in cur.fetchall()]
   
def get_object_contin_name(cur,obj):
    """
    Returns contin name of object number

    Parameters:
    -----------
    cur : MySQLdb cursor
    obj : str 
        object number
    """
    sql = ("select CON_AlternateName "
            "from contin "
            "join object on object.CON_Number = contin.CON_Number "
            "where object.OBJ_Name = %s " %(obj))
    try:
        cur.execute(sql)
        tmp = cur.fetchone()
        return tmp[0]
    except:
        return None

def get_object_contin(cur,obj):
    """
    Returns contin of object number

    Parameters:
    -----------
    cur : MySQLdb cursor
    obj : str 
        object number
    """
    sql = ("select contin.CON_Number "
            "from contin "
            "join object on object.CON_Number = contin.CON_Number "
            "where object.OBJ_Name = %s " %(obj))
    cur.execute(sql)
    return cur.fetchone()[0]


def get_contin_edges(cur,contin,start=0,end=None):
    """
    Returns list of edges for contin

    Parameters:
    -----------
    cur : MySQLdb cursor
    contin : str
        Contin number
    """
    if start or end:
        images = get_img_number(cur,start=start,end=end)
        images = ''.join(["'","','".join(images),"'"])
        sql = ("select OBJName1,OBJName2 from relationship "
                "join object as obj1 on obj1.OBJ_Name = relationship.OBJName1 "
                "join object as obj2 on obj2.OBJ_Name = relationship.OBJName2 "
                "where continNum = %s "
               "and (obj1.IMG_Number in (%s) "
               "or obj2.IMG_Number in (%s)) " %(contin,images,images))
    else:
        sql = ("select OBJName1,OBJName2 from relationship "
            "where continNum = %s" %contin)

    cur.execute(sql)
    return cur.fetchall()

def get_img_number(cur,**kwargs):
    """
    Returns list of image(sections) names from the NR an VC series

    Parameters:
    cur : MySQLdb cursors
    start : int
       First section number (default -1)
    end : int
       Last section number (default 1e6)
    """   
    args = aux.format.get_args(kwargs,{'start':-1,'end':1e6})
    sql = ("select IMG_Number "
           "from image "
           "where IMG_Series in ('NR','VC') " 
           "and IMG_SectionNumber >= %d "
           "and IMG_SectionNumber <= %d" %(int(args['start']),int(args['end'])))
    cur.execute(sql)
    return [i[0] for i in cur.fetchall()]

def convert_image_name(cur,imgNum):
    """
    Converts imgNum to image file name without extension, as used in adjacency2 table

    Parameters:
    -----------
    cur : MySQLdb cursors
    imgNum : str
        IMG_Number field entry for table image
    """
    sql = ("select IMG_File "
            "from image "
            "where IMG_Number = '%s' " %(imgNum))
    cur.execute(sql)
    return cur.fetchone()[0].split('.')[0]


def get_obj_area(cur):
    """
    Returns dictionary of (key,value) = (object_number,[image_number,area])

    Parameters:
    -----------
    cur : MySQLdb cursor

    """
    sql = ("select object,imgNum,area "
           "from dimensions ")
    cur.execute(sql)
    return dict([(a[0],{'image':a[1],'area':a[2]}) for a in cur.fetchall()])
           

def get_synapse_contins(cur,stype,**kwargs):
    """
    Returns list of synapse contin numbers

    Parameters:
    -----------
    cur : MySQLdb cursors
    stype : str
      synapses type; ['chemical','electrical']
    start : int
      beginning section number (default -1)
    end : int
      last section number (default 1e6)

    """
    args = aux.format.get_args(kwargs,{'start':None,'end':None,})
    if args['start'] or args['end']: 
        images = get_img_number(cur,start=args['start'],end=args['end'])
        images = ''.join(["'","','".join(images),"'"])
        sql = ("select synapsecombined.continNum "
               "from synapsecombined "
               "join object on object.CON_Number = synapsecombined.continNum "
               "where synapsecombined.type like '%s' "
               "and object.IMG_Number in (%s) " %(stype,images))
    else:
        sql = ("select synapsecombined.continNum "
               "from synapsecombined "
               "join object on object.CON_Number = synapsecombined.continNum "
               "where synapsecombined.type like '%s' " %(stype))			
    cur.execute(sql)
    return ','.join([str(a[0]) for a in cur.fetchall()])   

def get_presynapse_contins(cur,cell,start=0,end=None):
    """
    Returns lists of synapse contins wher cell is presynaptic

    Parameters:
    -----------
    cur : MySQLdb cursors
    cell: str
        Cell name
    start: int (optional)
        Start image index
    end: int (optional)
        End image index
    """
    if start or end:
        images = get_img_number(cur,start=start,end=end)
        images = ''.join(["'","','".join(images),"'"])
        sql = ("select synapsecombined.continNum "
               "from synapsecombined "
               "join object on object.CON_Number = synapsecombined.continNum "
               "where synapsecombined.type like 'chemical' "
               "and pre like '%%%s%%' "
               "and object.IMG_Number in (%s) " %(cell,images))
    else:
        sql = ("select synapsecombined.continNum "
               "from synapsecombined "
               "join object on object.CON_Number = synapsecombined.continNum "
               "where synapsecombined.type like 'chemical' " 
               "and pre like '%%%s%%' " %(cell))
        
    cur.execute(sql)
    return list(set([str(a[0]) for a in cur.fetchall()]))   

def get_postsynapse_contins(cur,cell,start=0,end=None):
    """
    Returns lists of synapse contins wher cell is postynaptic

    Parameters:
    -----------
    cur : MySQLdb cursors
    cell: str
        Cell name
    start: int (optional)
        Start image index
    end: int (optional)
        End image index
    """
    if start or end:
        images = get_img_number(cur,start=start,end=end)
        images = ''.join(["'","','".join(images),"'"])
        sql = ("select synapsecombined.continNum "
               "from synapsecombined "
               "join object on object.CON_Number = synapsecombined.continNum "
               "where synapsecombined.type like 'chemical' "
               "and post like '%%%s%%' "
               "and object.IMG_Number in (%s) " %(cell,images))
    else:
        sql = ("select synapsecombined.continNum "
               "from synapsecombined "
               "join object on object.CON_Number = synapsecombined.continNum "
               "where synapsecombined.type like 'chemical' " 
               "and post like '%%%s%%' " %(cell))
        
    cur.execute(sql)
    return list(set([str(a[0]) for a in cur.fetchall()]))   

def get_gapjunction_contins(cur,cell,start=0,end=None):
    """
    Returns lists of synapse contins wher cell is a gap junction

    Parameters:
    -----------
    cur : MySQLdb cursors
    cell: str
        Cell name
    start: int (optional)
        Start image index
    end: int (optional)
        End image index
    """
    if start or end:
        images = get_img_number(cur,start=start,end=end)
        images = ''.join(["'","','".join(images),"'"])
        sql = ("select synapsecombined.continNum "
               "from synapsecombined "
               "join object on object.CON_Number = synapsecombined.continNum "
               "where synapsecombined.type like 'electrical' "
               "and (post like '%%%s%%' or pre like '%%%s%%') "
               "and object.IMG_Number in (%s) " %(cell,cell,images))
    else:
        sql = ("select synapsecombined.continNum "
               "from synapsecombined "
               "join object on object.CON_Number = synapsecombined.continNum "
               "where synapsecombined.type like 'electrical' " 
               "and (post like '%%%s%%' or pre like '%%%s%%' " %(cell,cell))
        
    cur.execute(sql)
    return list(set([str(a[0]) for a in cur.fetchall()]))   

def get_synapse_data_by_contin(cur,contin):
    """
    Returns synapse data for given contin
    Each row is a single section of the synapse
    Row format: [section_number,preobj,[post_obj]]

    Parameters:
    -----------
    cur : MySQLdb cursor
    contin : str
        Contin number
    """
    sql = ("select IMG_Number,fromObj,toObj "
            "from object "
            "where CON_Number = %s " %(contin))
    cur.execute(sql)
    try:
        return [(a[0],a[1],a[2].split(',')) for a in cur.fetchall()]
    except:
        return None

def get_synapse_data(cur,stype,**kwargs):
    """
    Returns list of synapses
    [pre_cell,post_cell,number_of_sections,continNum,series]
    
    Parameters:
    -----------
    cur : MySQLdb cursor
    stype : str
      synapse type; ['chemical','electrical']
    start : int
      beginning section number (default -1)
    end : int
      last section number (default 1e6)    
    
    """
    contins = get_synapse_contins(cur,stype,**kwargs)
    sql = ("select pre,post,sections,continNum,series "
           "from synapsecombined "
           "where continNum in (%s) " %(contins))	
    cur.execute(sql)
    return cur.fetchall()
    
def get_adjacency_data(cur):
    """
    Return list of adjacency data,
    [cell1,cell2,amount_of_contact,image_number]

    Parameters:
    -----------
    cur : MySQLdb cursor

    """
    sql = ("select pre,post,weight,imgNum "
           "from adjacency2")
    cur.execute(sql)
    return [(a[0],a[1],int(a[2]),a[3]) for a in cur.fetchall()]

def get_adjacency_from_layer(cur,layer):
    """
    Return list of adjacency data from given layer
    [cell1,cell2,amount_of_contact,image_number]

    Parameters:
    -----------
    cur : MySQLdb cursor
    
    """
    sql = ("select pre,post,weight "
           "from adjacency2 "
           "where imgNum = '%s'" %layer)
    cur.execute(sql)
    return [(a[0],a[1],int(a[2])) for a in cur.fetchall()]



def get_adjacency_cells(cur):
    """
    Returns list of cells in the adjacency data

    Parameters:
    -----------
    cur : MySQLdb cursor

    """
    sql = ("select distinct(pre) from adjacency2 "
           "union "
           "select distinct(post) from adjacency2")
    cur.execute(sql)
    return [a[0] for a in cur.fetchall()]


def get_touch_density(cur,key="pixels_norm"):
    """
    Returns list of 'touch_density' measures;
    [cell1,cell2,touch_density]

    Parameters:
    -----------
    cur : MySQLdb cursor
    key : str
     touch density measure [pixels,pixels_norm,segments,segments_norm]
    
    """
    sql = ("select pre,post,%s "
           "from touch_density" %key)
    cur.execute(sql)
    return [(a[0],a[1],float(a[2])) for a in cur.fetchall()]

def maxAnterior(cur,cell,contin=None):
    """
    Returns max anterior section of cell

    Parameters:
    -----------
    cur : MySQLdb cursor
    cell : str
      Name of cell
    contin : str (optional, defautl = None)
      Contin number
    """
    if contin:
        sql = ("select min(image.IMG_SectionNumber) "
               "from image "
               "join object on image.IMG_Number = object.IMG_Number "
               "join contin on contin.CON_Number=object.CON_Number "
               "where (contin.CON_AlternateName like '%%%s%%' ) "
               "and (contin.CON_Number = %d) "
               #%(cell,contin))
               %(cell,int(contin)))
    else:
        sql = ("select min(image.IMG_SectionNumber) "
               "from image "
               "join object on image.IMG_Number = object.IMG_Number "
               "join contin on contin.CON_Number=object.CON_Number "
               "where contin.CON_AlternateName like '%%%s%%' "
               %cell)    
    cur.execute(sql)
    return cur.fetchone()[0]

def maxPosterior(cur,cell,contin=None):
    """
    Returns the max posterior sections of cell

    Parameters:
    -----------
    cur : MySQLdb cursor
    cell : str
      name of cell
    contin : str (optional, defautl = None)
      Contin number
    """

    if contin:
        sql = ("select max(image.IMG_SectionNumber) "
               "from image "
               "join object on image.IMG_Number = object.IMG_Number "
               "join contin on contin.CON_Number=object.CON_Number "
               "where (contin.CON_AlternateName like '%%%s%%' ) "
               "and (contin.CON_Number = %d) "
               %(cell,int(contin)))
    else:
        sql = ("select max(image.IMG_SectionNumber) "
               "from image "
               "join object on image.IMG_Number = object.IMG_Number "
               "join contin on contin.CON_Number=object.CON_Number "
               "where contin.CON_AlternateName like '%%%s%%' "
               %cell)

    cur.execute(sql)
    return cur.fetchone()[0]

def neuron_cylinder(cur,name):
    """
    Returns the cylindrical locations of neuron segment centroids
    [cell_name,radius,phi,z]

    Parameters:
    ----------
    cur : MySQLdb
    name : str
        Cell name
    """
    sql = ("select radialPharynx.distance,radialPharynx.phi,image.IMG_SectionNumber "
            "from radialPharynx "
            "join object on object.OBJ_Name=radialPharynx.OBJ_Name "
            "join contin on contin.CON_Number=object.CON_Number "
            "join image on image.IMG_Number=object.IMG_Number "
            "where contin.CON_AlternateName = '%s'" %name)

    cur.execute(sql)
    data = [[int(a[0]),float(a[1]),int(a[2])] for a in cur.fetchall()]
    return data     
    
def gap_cylinder(cur,dr=0,dphi=0.1):
    """
    Return the cylindrical coordinates of the gap junction
    [pre_cell,post_cell,radius,phi,z]

    Parameters:
    -----------
    cur : MySQLdb cursor
    dr : int
      shift for radius values (default 0)
    dphi : float
      shift for phi values (default 0.1)

    """
    sql = ("select synapsecombined.pre,"
           "synapsecombined.post,"
           "radialPharynx.distance,"
           "radialPharynx.phi,"
           "image.IMG_SectionNumber "
           "from synapsecombined "
           "join radialPharynx on radialPharynx.OBJ_Name = synapsecombined.mid "
           "join object on object.OBJ_Name = radialPharynx.OBJ_Name "
           "join image on image.IMG_Number=object.IMG_Number "
           "where synapsecombined.type='electrical'"
           )

    cur.execute(sql)
    gap = []
    for a in cur.fetchall():
        r = int(a[2])
        phi = float(a[3])
        z = int(a[4])
        gap.append([a[0],a[1],r,phi,z])
        
    return gap    

def syn_cylinder(cur):
    """
    Return the cylindrical coordinates of synapses
    [pre_cell,post_cell,radius,phi,z]

    Parameters:
    -----------
    cur : MySQLdb cursor

    """
    sql = ("select synapsecombined.pre,"
           "synapsecombined.post,"
           "radialPharynx.distance,"
           "radialPharynx.phi,"
           "image.IMG_SectionNumber "
           "from synapsecombined "
           "join radialPharynx on radialPharynx.OBJ_Name = synapsecombined.mid "
           "join object on object.OBJ_Name = radialPharynx.OBJ_Name "
           "join image on image.IMG_Number=object.IMG_Number "
           "where synapsecombined.type='chemical'"
           )

    cur.execute(sql)
    syn = []
    for a in cur.fetchall():
        r = int(a[2])
        phi = float(a[3])
        z = int(a[4])
        pre = a[0]
        for post in a[1].split(','):
            syn.append([pre,post,r,phi,z])
    return syn

def get_cell_genes(cur,cell):
    """
    Return the genes expressed in cell

    Parameters
    ----------
    cur : MySQLdb cursor
    cell : str
      Name of cell
    """
    sql = ("select genename from WB.association "
           "join WB.cells on "
           "WB.association.idcell = WB.cells.idcells "
           "where WB.cells.name like '%s'"
           %cell)
    cur.execute(sql)
    return [a[0] for a in cur.fetchall()]

def get_object_xy_in_layer(cur,cell,layer):
    """
    Returns (object_num,x,y) coordinates of cell object in layer

    Parameters
    ----------
    cur : MySQLdb cursor
    cell : str
      Name of cell
    layer : str
      Layer/image name
    """

    sql = ("select OBJ_Name,OBJ_X,OBJ_Y "
           "from object "
           "join contin on "
           "contin.CON_Number = object.CON_Number "
           "where contin.CON_AlternateName like '%s' "
           "and object.IMG_Number like '%s'"
           %(cell,layer))
    cur.execute(sql)
    return [(int(a[0]),int(a[1]),int(a[2])) for a in cur.fetchall()]
            
def get_object_xyz(cur,obj):
    """
    Returns (x,y,z) coordinates of object

    Parameters:
    -----------
    cur : MySQLdb cursor
    obj : str
        object name
    """
    sql = ("select OBJ_X,OBJ_Y,image.IMG_SectionNumber "
            "from object "
            "join image on image.IMG_Number = object.IMG_Number "
            "where OBJ_Name = '%s' " %(obj))
    cur.execute(sql)
    return cur.fetchone()

def get_synapse_from_layer(cur,layer):
    """
    Returns synapses in given layer list of tuples:
    (mid_object,x-pos,y-pos)
    
    Parameters
    ----------
    cur : MySQLdb cursor
    layer : str
      Name of layer
    
    """

    sql = ("select synapsecombined.mid,"
           "object.OBJ_X,"
           "object.OBJ_Y "
           #"image.IMG_SectionNumber "
           "from synapsecombined "
           "join object on object.OBJ_Name=synapsecombined.mid "
           "join image on image.IMG_Number=object.IMG_Number "
           "where object.IMG_Number='%s' "
           #"and synapsecombined.type='chemical'"
           %layer)
    cur.execute(sql)
    return [(int(a[0]),int(a[1]),int(a[2])) for a in cur.fetchall()]

def get_objects_in_layer(cur,layer):
    """
    Return object coordinates in layer
    """
    sql = ("select preObj from adjacency2 where imgNum = '%s' "
            "union "
            "select postObj from adjacency2 where imgNum = '%s'"
            %(layer,layer))
    
    cur.execute(sql)
    return list(set([a[0] for a in cur.fetchall()]))

def get_object_adjacency(cur,obj):
    """
    Returns the adjacency of the object number

    Parameters
    ----------
    cur : MySQLdb cursor
    obj : int
       Object number

    """
    sql = ("select pre from adjacency2 "
           "where postObj = %s "
           "union "
           "select post from adjacency2 "
           "where preObj = %s "
           %(obj,obj)
           )
    try:
        cur.execute(sql)
        return [a[0] for a in cur.fetchall()]
    except:
        return None

def order_synapses_by_section_number(cur,stype='chemical'):
    """
    Returns synapses ordered by section number. Row format
    [pre,post,sections,mid obj,continNum,IMG_Number,IMG_SectionNumber,X,Y]
    
    Parameters
    ----------
    cur : MySQLdb cursor
    stype : str
      synapse type : 'chemical' or 'electrical'

    """

    sql = ("select pre,post,sections,mid,continNum,"
           "object.IMG_Number,image.IMG_SectionNumber, "
           "object.OBJ_X,object.OBJ_Y "
           "from synapsecombined "
           "join object on synapsecombined.mid = object.OBJ_Name "
           "join image on image.IMG_Number = object.IMG_Number "
           "where synapsecombined.type = '%s' "
           "order by image.IMG_SectionNumber"
           %stype)
    cur.execute(sql)
    return cur.fetchall()

def get_contin_xyz(cur,contin):
    """"
    Returns objects and associated coordinates for contin 
    [OBJ_Number,OBJ_X,OBJ_Y,IMG_SectionNumber,IMG_File]

    Parameters
    ----------
    cur : MySQLdb cursor
    contin : int 
        contin number

    """

    sql = ("select OBJ_Name,OBJ_X,OBJ_Y,IMG_SectionNumber,image.IMG_File " 
            "from object "
            "join image on object.IMG_Number = image.IMG_Number "
            "where CON_Number = %d"
            %contin)
    cur.execute(sql)
    return cur.fetchall()

