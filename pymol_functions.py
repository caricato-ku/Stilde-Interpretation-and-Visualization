from pymol import cmd,preset,util
import pandas as pd
import numpy as np
from pymol.cgo import *

#Keep track of times elec_mag is called
count=1

def loadCSV(filename):
    data = pd.read_csv(filename)
    if data.empty:
        print("Dataframe is empty")
    else:
        print(data[['S','Electric Magnitude', 'Magnetic Magnitude','Cosine of Angle']])
    return data
cmd.extend("loadCSV", loadCSV)

def newLoad(filename):
    cmd.load(filename)
    # Everything ball and stick
    preset.ball_and_stick(selection="all", mode=1)

    # Change carbon color
    cmd.color("gray",selection="elem c")
    return
cmd.extend("newLoad",newLoad)

def str_to_list(string,internalType="float"):
   """PyMol reads everything as strings, need to convert
      convert to list for Python
   """
   if type(string) is list:
       return string
   newList=string.strip('][').split(',')
   if internalType=="float":
       newList=[float(n) for n in newList]
   return newList


def elec_mag(elec_end,mag_end,elec_scaling_factor=7, mag_scaling_factor=7, 
             elec_start=[0.0,0.0,0.0],mag_start=[0.0,0.0,0.0]):
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
    #Just calls the other function with appropriate args,
    #avoids repeating function body
    elec_mag(elec_end,mag_end,elec_start=elec_start,mag_start=mag_start, elec_scaling_factor=elec_scale, mag_scaling_factor=mag_scale)
    return
cmd.extend("elec_mag_fromAtom", elec_mag_fromAtom)


def select_vectors(index, df):
    index = int(index)
    cart=['X','Y','Z']
    ##Uses list comprehension to express more succinctly
    elecVec = [df.iloc[index]['Electric'+x] for x in cart]
    magVec = [df.iloc[index]['Magnetic'+x] for x in cart]
    vecList = [elecVec, magVec]
    return vecList
cmd.extend("select_vectors",select_vectors)

def multiple_vectors(indices, df, fromAtom=False):
    for index in indices:
        vec = select_vectors(index, df)
        if fromAtom == True:
            elec_mag_fromAtom(vec[0], vec[1])
        else:
            elec_mag(vec[0], vec[1])

cmd.extend("multiple_vectors", multiple_vectors)

def createSphere(pos, radius=10.0, color = 'Yellow',transparency=.5):
    cmd.set("cgo_sphere_quality", 4)
    radius=float(radius)
    pos = str_to_list(pos)
    obj = [cgo.SPHERE] + pos + [radius]
    cmd.load_cgo(obj,'s1',0)
    cmd.color(color,selection='s1')
    cmd.set("cgo_transparency",value=transparency,selection="s1")
cmd.extend("createSphere",createSphere)
