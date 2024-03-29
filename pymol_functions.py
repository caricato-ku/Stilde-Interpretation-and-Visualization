from pymol import cmd,preset,util,cgo
import pandas as pd
import numpy as np
from copy import deepcopy

#Keep track of times elec_mag is called
count=1

def loadCSV(filename):
    """
    Loads and display a CSV file via Panda dataframe

    :param filename: Name of the file to be opened
    :type filename: string


    :return: Dataframe of the CSV opened
    :rtype: Dataframe
   """
    data = pd.read_csv(filename)
    if data.empty:
        print("Dataframe is empty")
    else:
        print(data[['S','Electric Magnitude', 'Magnetic Magnitude','Cosine of Angle']])
    return data
cmd.extend("loadCSV", loadCSV)


def newLoad(filename):

    """
    Loads a new molecule file into PyMol, setting carbons to gray and using ball and stick representation

    :param filename: Name of molecule file to be opened

    :return: None
    """
    cmd.load(filename)
    # Everything ball and stick
    preset.ball_and_stick(selection="all", mode=1)

    # Change carbon color
    cmd.color("gray",selection="elem c")
    return
cmd.extend("newLoad",newLoad)


def str_to_list(string,internalType="float"):
    """
        :param string: Object to be converted to list of floats
        :type string: string, float, or np.ndarray

        :param internalType: String indicating the type of list to be converted to
        :type internalType: string (,optional -) defaults to float

        :return: Converted list of objectss
        :rtype: List of floats
    """
    if type(string) in (list,np.ndarray):
       return list(string)
    newList=string.strip('][').split(',')
    if internalType=="float":
       newList=[float(n) for n in newList]
    return newList

def elec_mag(elec_end,mag_end,elec_scale=7, mag_scale=7, 
             elec_start=[0.0,0.0,0.0],mag_start=[0.0,0.0,0.0],
             use_lab=True):
    """
    Main driver function for drawing arrows in PyMol. Calls cgo_arrow function with correct parameters

    :param elec_end: Endpoint of electric vector
    :type elec_end: List of floats

    :param mag_end: Endpoint of magnetic vector
    :type mag_end: List of floats

    :param elec_scale: Scaling factor for electric vector
    :type elec_scale: int, optional - defaults to 7

    :param mag_scale: Scaling factor for magnetic vector
    :type mag_scale: int, optional - defaults to 7

    :param elec_start: Starting point for electric vector - defaults to coordinate origin [0,0,0]
    :type elec_start: List of floats

    :param mag_start: Starting point for magnetic vector - defaults to coordinate origin [0,0,0]
    :type mag_start: List of floats


    :return: None
        """
    global count
    elec_end=str_to_list(elec_end)
    temp_elec_end = elec_end.copy()
    if elec_start!='sele':
        elec_start=str_to_list(elec_start)
    #Use count to make a unique name for each arrow object,
    #allows multiple elec_mags to be spawned
    cgo_arrow(origin=elec_start,endpoint=temp_elec_end,type="electric",name=str(count), 
              scaling=elec_scale,use_lab=use_lab)

    mag_end=str_to_list(mag_end)
    temp_mag_end = mag_end.copy()
    if mag_start!='sele':
        mag_start=str_to_list(mag_start)

    
    cgo_arrow(origin=mag_start,endpoint=temp_mag_end,type="magnetic",name=str(count),
              scaling=mag_scale,use_lab=use_lab)
    cmd.group(f"stilde{count}",members=f"electric{count} magnetic{count}")

    count+=1
    return
cmd.extend("elec_mag",elec_mag)


def elec_mag_fromAtom(elec_end,mag_end, elec_scale=7, mag_scale=7, elec_start='sele',mag_start='sele'):
    """
    Similar functionality to elec_mag, but draws arrow using atom as origin point

    :param elec_end: Endpoint of electric vector
    :type elec_end: List of floats

    :param mag_end: Endpoint of magnetic vector
    :type mag_end: List of floats

    :param elec_scale: Scaling factor for electric vector
    :type elec_scale: int, optional - defaults to 7

    :param mag_scaling_factor: Scaling factor for magnetic vector
    :type mag_scaling_factor: int, optional - defaults to 7

    :param elec_start: Starting point for electric vector - defaults to selected atom xyz coordinates
    :type elec_start: List of floats

    :param mag_start: Starting point for magnetic vector - defaults to selected atom xyz coordinates
    :type mag_start: List of floats


    :return: None
    """
    #Just calls the other function with appropriate args,
    #avoids repeating function body
    elec_mag(elec_end,mag_end,elec_start=elec_start,mag_start=mag_start, elec_scale=elec_scale, mag_scale=mag_scale)
    return
cmd.extend("elec_mag_fromAtom", elec_mag_fromAtom)


def select_vectors(index, df, fromAtom=False):
    """
    Pulls vector data from dataframe based on a given index. Automatically calls elec_mag to draw arrows

    :param index: Index of dataframe where vector data will be selected from
    :type index: Int

    :param df: Dataframe containing vector data loaded previously
    :type df: Dataframe

    :param fromAtom: Boolean indicating whether or not the vector will be drawn using an atom as origin, defaults to False
    :type fromAtom: Boolean

    :return: List of lists containing electric vector at index 0, magnetic vector at index 1
    :rtype: List of ints
    """
    index = int(index)
    cart=['X','Y','Z']
    ##Uses list comprehension to express more succinctly
    elecVec = [df.iloc[index]['Electric'+x] for x in cart]
    magVec = [df.iloc[index]['Magnetic'+x] for x in cart]
    vecList = [elecVec, magVec]
    if fromAtom is False:
        elec_mag(elecVec, magVec,elec_scale=2,mag_scale=2)
    else:
        elec_mag_fromAtom(elecVec,magVec,elec_scale=2,mag_scale=2)
    return vecList
cmd.extend("select_vectors",select_vectors)

def multiple_vectors(indices, df, fromAtom=False):

    """
    Similar functionality to select_vectors, but calls select_vectors in a loop
    for a given list of indices to draw multiple arrows

    :param indices: Indices of dataframe where vector data will be selected from
    :type index: List of ints

    :param df: Dataframe containing vector data loaded previously
    :type df: Dataframe

    :param fromAtom: Boolean indicating whether or not the vector will be drawn using an atom as origin, defaults to False
    :type fromAtom: Boolean
    :return: None
    """

    for index in indices:
        vec = select_vectors(index, df, fromAtom)

cmd.extend("multiple_vectors", multiple_vectors)

def createSphere(pos, radius=1.0, color = 'Yellow',transparency=.5):

    """
    Draws a sphere representative of s-tilde magnitude

    :param pos: Indicates the x,y,z coordinates for the sphere to be drawn at
    :type pos: List of floats

    :param radius: Radius of sphere - default 1.0
    :type radius: Float

    :param color: Color of sphere - default yellow
    :type color: String

    :param transparency: transparency value of sphere - default .5
    :type transparency: Float

    Return value: None
    """
    cmd.set("cgo_sphere_quality", 4)
    radius=float(radius)
    pos = str_to_list(pos)
    obj = [cgo.SPHERE] + pos + [radius]
    cmd.load_cgo(obj,'s1',0)
    cmd.color(color,selection='s1')
    cmd.set("cgo_transparency",value=transparency,selection="s1")
cmd.extend("createSphere",createSphere)

def checkVecs(pairs,gauge="VE"):
    """Use gaugeComp to check given pairs

    """
    vectors0=pd.read_csv(f"vector_{gauge}0N",delim_whitespace=True, header=None)    
    vectors0.set_index([0,1],inplace=True)    
    vectors1=pd.read_csv(f"vector_{gauge}1N",delim_whitespace=True, header=None)    
    vectors1.set_index([0,1],inplace=True)

    for i,a in pairs:
        vecs0=np.array(vectors0.loc[i,a])
        elec0,mag0=vecs0[1:4],vecs0[4:]
        vecs1=np.array(vectors1.loc[i,a])
        elec1,mag1=vecs1[1:4],vecs1[4:]
        if "E" in gauge:
            elec1-=elec0
            elec1*=-1.0
        elif "M" in gauge:
            mag1-=mag0
            mag1*=-1.0
        
        result={gauge:[elec1,mag1]}
        gaugeComp(result)
cmd.extend("checkVecs",checkVecs)


def gaugeComp(vecs,freq=0.077357,max_len=4.0,use_lab=False):
    """Compares stilde vectors between gauges

    Receives a dict with lg,mvg-m, and mvg-e lists of vectors.
    Draws them after rescaling relative to the longest. MVG vectors each 
    receive a factor of (1/freq)^(0.5). 
    """
    vecs=deepcopy(vecs)
    print(vecs)
    inv_freq=(1.0/freq)**.5
    for gauge in vecs:
        if gauge in ('VM','VE'):
            vecs[gauge][0]=[v*inv_freq for v in vecs[gauge][0]]
            vecs[gauge][1]=[v*inv_freq for v in vecs[gauge][1]]

    ##Determine longest elec/magnetic, scale relative to these
    lengths_E={gauge:sum([v**2 for v in vecs[gauge][0]])**.5 for gauge in vecs}
    lengths_M={gauge:sum([v**2 for v in vecs[gauge][1]])**.5 for gauge in vecs}

    stildes={gauge:np.dot(vecs[gauge][0],vecs[gauge][1]) for gauge in vecs}
    ang=np.degrees([np.arccos(stildes[gauge]/(lengths_E[gauge]*lengths_M[gauge])) for gauge in vecs])

    max_E=max(lengths_E.values())
    max_M=max(lengths_M.values())    

    for gauge in vecs:
        Escale=(max_len/max_E)
        Mscale=(max_len/max_M)

        elec_mag(vecs[gauge][0],vecs[gauge][1],elec_scale=Escale,mag_scale=Mscale,use_lab=use_lab)
    print(f"Stilde={stildes}\n")
    print(f"Electric={lengths_E}\n")
    print(f"Magnetic={lengths_M}\n")
    print(f"Angles={ang}\n")
    return
cmd.extend("gaugeComp",gaugeComp)

def polar_plot(tensor=None, Ntheta=150, Nphi=150, origin=None,
               scale=1.0, pos_color="Blue", neg_color="Orange",style="loop"):
    """Generate a polar plot from the passed in tensor

    Ntheta/Nphi give how many increments of theta/phi to sample.
    origin lets you select a displacement of the polar plot.
    scale is a factor applied to r for all points

    Can also select colors for positive/negative portion 
    (current default matches Autschbach) as well as style of the plot
    ('points' just plots vertices, 
     'lines' draws small segments between pairs of points,
     'loops' draws a continuous line from the first to last point)
    """
    if tensor is None:
        tensor = np.array([1.0,1.0,1.0,0.0,0.0,0.0])
    if origin is None:
        origin = [0.0,0.0,0.0]

    pos=[]
    neg=[]

    theta_increment=(2*np.pi)/Ntheta
    phi_increment=(np.pi)/Nphi

    # For 'loops' and 'lines' style, alternating
    # between forward and reverse loops over phi
    # avoids the next linked vertex being far away.
    forward = range(Nphi)
    reverse = range(Nphi,0,-1)

    for i in range(Ntheta): 
        theta=i*theta_increment
        if i%2==0:
            phi_range=forward
        else:
            phi_range=reverse
        for j in phi_range:
            phi=j*phi_increment
            r=scale*calc_r(theta,phi,tensor)
            if r>=0.0:
                pos.append((r,theta,phi))
            else:
                neg.append((np.abs(r),theta,phi))
    
    cart_pos = convert_cartesian(pos)        
    cart_neg = convert_cartesian(neg)        

    if style.lower() == "loop":
        begin = [cgo.BEGIN, cgo.LINE_LOOP]
    elif style.lower() == "points":
        begin = [cgo.BEGIN, cgo.POINTS]
    elif style.lower() == "lines":
        begin = [cgo.BEGIN, cgo.LINES]
    else:
        print("Unknown style, defaulting to 'loop'")
        begin = [cgo.BEGIN, cgo.LINE_LOOP]

    # Write out list of vertices to make lines. For all but the first
    # and last vertex need to repeat (e.g. v1->v2, v2->v, ...)
    obj_pos = [[cgo.VERTEX,x[0]+origin[0],x[1]+origin[1],x[2]+origin[2]] for x in cart_pos]
    #flatten the nested list
    obj_pos = [vert for verts in obj_pos for vert in verts]
    #Append group specifiers for CGO object
    obj_pos = begin + obj_pos + [cgo.END]

    obj_neg = [[cgo.VERTEX,x[0],x[1],x[2]] for x in cart_neg]
    obj_neg = [vert for verts in obj_neg for vert in verts]
    obj_neg = begin + obj_neg + [cgo.END]

    cmd.load_cgo(obj_pos,"positive_polar",0)
    cmd.color(pos_color,selection='positive_polar')
    cmd.load_cgo(obj_neg,"negative_polar",0)
    cmd.color(neg_color,selection='negative_polar')
    return
cmd.extend("polar_plot",polar_plot)

def calc_r(theta=0.0,phi=0.0,tensor=None):
    """Evaluate r from the polar function"""
    sp=np.sin(phi)
    st=np.sin(theta)
    cp=np.cos(phi)
    ct=np.cos(theta)

    r=(tensor[0]*sp**2*ct**2
      +tensor[1]*sp**2*st**2
      +tensor[2]*cp**2
      +2*tensor[3]*st*ct*sp**2
      +2*tensor[4]*ct*sp*cp
      +2*tensor[5]*st*sp*cp)

    return r
     
def convert_cartesian(polar_coords):
    """Convert polar coordinates to Cartesian

    Expects `polar_coords` to be a list of (r,theta,phi) tuples.
    """
    return [(r*np.cos(theta)*np.sin(phi), r*np.sin(theta)*np.sin(phi), r*np.cos(phi))
            for r,theta,phi in polar_coords]


def coord_axes(x_scale=1.0,y_scale=1.0,z_scale=1.0, loc=None):
    """Display individually scalable coordinate axes

    Adapted from https://pymolwiki.org/index.php/Axes
    """
    if loc==None:
        loc = [1.0,1.0,1.0]

    w = 0.06 # cylinder width 
    l = 0.75 # cylinder length
    h = 0.25 # cone height
    d = w * 1.618 # cone base diameter
    
    x=x_scale*l
    y=y_scale*l
    z=z_scale*l

    x_cone_start = loc.copy()
    x_cone_start[0] += x
    x_cone_end = x_cone_start.copy()
    x_cone_end[0] += h

    y_cone_start = loc.copy()
    y_cone_start[1] += y
    y_cone_end = y_cone_start.copy()
    y_cone_end[1] += h

    z_cone_start = loc.copy()
    z_cone_start[2] += z
    z_cone_end = z_cone_start.copy()
    z_cone_end[2] += h

    red = [1.0,0.0,0.0]
    green = [0.0,1.0,0.0]
    blue = [0.0,0.0,1.0]

    #CYLINDER is defined by a top center (x1,y1,z1), bottom center (x2,y2,z2), radius (r),
    #start color (r1,g1,b1), and end color (r2,g2,b2) (allows for a color gradient)

    #CONE is top center (x1,y1,z1), bottom center (x2,y2,z2), top radius (r1), bottom radius (r2),
    #start color (r1,g1,b1), end color (r2,g2,b2), and 0/1 for open/closed top/bottom       

    x_arrow = ([CYLINDER] + loc + x_cone_start + [w] + 2*red
               +[CONE] + x_cone_start + x_cone_end + [d] + [0.0] + 2*red + [1.0, 1.0]) 
    y_arrow = ([CYLINDER] + loc + y_cone_start + [w] + 2*green
               +[CONE] + y_cone_start + y_cone_end + [d] + [0.0] + 2*green + [1.0, 1.0]) 
    z_arrow = ([CYLINDER] + loc + z_cone_start + [w] + 2*blue
               +[CONE] + z_cone_start + z_cone_end + [d] + [0.0] + 2*blue + [1.0, 1.0])
    obj = x_arrow + y_arrow + z_arrow 

    cmd.load_cgo(obj, 'axes')
    return
cmd.extend("coord_axes",coord_axes)
