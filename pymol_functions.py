from pymol import cmd,preset,util
import pandas as pd
import numpy as np
from pymol.cgo import *

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
        :param string: String object to be converted to list of floats
        :type string: string

        :param internalType: String indicating the type of list to be converted to
        :type internalType: string (,optional -) defaults to float

        :return: Converted list of objectss
        :rtype: List of floats
    """
    if type(string) is list:
       return string
    newList=string.strip('][').split(',')
    if internalType=="float":
       newList=[float(n) for n in newList]
    return newList

def elec_mag(elec_end,mag_end,elec_scaling_factor=7, mag_scaling_factor=7, elec_start=[0.0,0.0,0.0],mag_start=[0.0,0.0,0.0]):
    """
    Main driver function for drawing arrows in PyMol. Calls cgo_arrow function with correct parameters

    :param elec_end: Endpoint of electric vector
    :type elec_end: List of floats

    :param mag_end: Endpoint of magnetic vector
    :type mag_end: List of floats

    :param elec_scaling_factor: Scaling factor for electric vector
    :type elec_scaling_factor: int, optional - defaults to 7

    :param mag_scaling_factor: Scaling factor for magnetic vector
    :type mag_scaling_factor: int, optional - defaults to 7

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
    cgo_arrow(origin=elec_start,endpoint=temp_elec_end,type="electric",name=str(count), scaling=elec_scaling_factor)

    mag_end=str_to_list(mag_end)
    temp_mag_end = mag_end.copy()
    if mag_start!='sele':
        mag_start=str_to_list(mag_start)

    cgo_arrow(origin=mag_start,endpoint=temp_mag_end,type="magnetic",name=str(count),scaling=mag_scaling_factor)
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

    :param elec_scaling_factor: Scaling factor for electric vector
    :type elec_scaling_factor: int, optional - defaults to 7

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
    elec_mag(elec_end,mag_end,elec_start=elec_start,mag_start=mag_start, elec_scaling_factor=elec_scale, mag_scaling_factor=mag_scale)
    return
cmd.extend("elec_mag_fromAtom", elec_mag_fromAtom)


#Function can't (or at least shouldn't) access variables from inside
#other functions.
##If you want to make this default, can set df=None and call
##loadCSV with a set file in this case.
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
        elec_mag(elecVec, magVec)
    else:
        elec_mag_fromAtom(elecVec,magVec)
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
