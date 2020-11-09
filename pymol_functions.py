from pymol import cmd,preset,util
import pandas as pd


def loadCSV(filename):
    data = pd.read_csv(filename)
    if data.empty:
        print("Dataframe is empty")
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


def elec_mag(elec_end,mag_end,elec_start=[0.0,0.0,0.0],mag_start=[0.0,0.0,0.0]):

    tempEndElec = elec_end.copy()
    elec_end=str_to_list(tempEndElec)
    elec_start=str_to_list(elec_start)
    cgo_arrow(origin=elec_start,endpoint=elec_end,type="electric")

    tempEndMag = mag_end.copy()
    mag_end=str_to_list(tempEndMag)
    mag_start=str_to_list(mag_start)
    cgo_arrow(origin=mag_start,endpoint=mag_end,type="magnetic")

    return
cmd.extend("elec_mag",elec_mag)


def elec_mag_fromAtom(elec_end,mag_end,elec_start='sele',mag_start='sele'):
    tempEndElec = elec_end.copy()
    elec_end=str_to_list(tempEndElec)
    cgo_arrow(origin=elec_start,endpoint=elec_end,type="electric")

    tempEndMag = mag_end.copy()
    mag_end=str_to_list(tempEndMag)
    cgo_arrow(origin=mag_start,endpoint=mag_end,type="magnetic")
    return
cmd.extend("elec_mag_fromAtom", elec_mag_fromAtom)


def select_vectors(index, df = data):
    index = int(index)
    elecVec = [df.iloc[index]['ElectricX'], df.iloc[index]['ElectricY'], df.iloc[index]['ElectricZ']]
    magVec = [df.iloc[index]['MagneticX'], df.iloc[index]['MagneticY'], df.iloc[index]['MagneticZ']]
    vecList = [elecVec, magVec]
    return vecList
cmd.extend("select_vectors",select_vectors)
