from pymol import cmd,preset,util
import pandas as pd

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


def elec_mag(elec_end,mag_end,elec_start=[0.0,0.0,0.0],mag_start=[0.0,0.0,0.0]):
    elec_end=str_to_list(elec_end)    
    elec_start=str_to_list(elec_start)
    cgo_arrow(atom1=elec_start,endpoint=elec_end,type="electric")

    mag_end=str_to_list(mag_end)
    mag_start=str_to_list(mag_start)
    cgo_arrow(atom1=mag_start,endpoint=mag_end,type="magnetic")
    return
cmd.extend("elec_mag",elec_mag)

def select_vectors(filename, index):
    data = pd.read_csv(filename)
    elecVec = [data.iloc[index]['ElectricX'], data.iloc[index]['ElectricY'], data.iloc[index]['ElectricZ']]
    magVec = [data.iloc[index]['MagneticX'], data.iloc[index]['MagneticY'], data.iloc[index]['MagneticZ']]
    vecList = [elecVec, magVec]
    return vecList
cmd.extend("select_vectors",select_vectors)
