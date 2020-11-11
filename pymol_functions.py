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
    elec_end=str_to_list(elec_end)    
    if elec_start!='sele':
        elec_start=str_to_list(elec_start)
    cgo_arrow(atom1=elec_start,endpoint=elec_end,type="electric")

    mag_end=str_to_list(mag_end)
    if mag_start!='sele':
        mag_start=str_to_list(mag_start)
    cgo_arrow(atom1=mag_start,endpoint=mag_end,type="magnetic")
    return
cmd.extend("elec_mag",elec_mag)


def elec_mag_fromAtom(elec_end,mag_end,elec_start='sele',mag_start='sele'):
    #Just calls the other function with appropriate args,
    #avoids repeating function body
    elec_mag(elec_end,mag_end,elec_start=elec_start,mag_start=mag_start)
    return
cmd.extend("elec_mag_fromAtom", elec_mag_fromAtom)


#Function can't (or at least shouldn't) access variables from inside
#other functions.
##If you want to make this default, can set df=None and call 
##loadCSV with a set file in this case. 
def select_vectors(index, df):
    index = int(index)
    cart=['X','Y','Z']
    ##Uses list comprehension to express more succinctly
    elecVec = [df.iloc[index]['Electric'+x] for x in cart]
    magVec = [df.iloc[index]['Magnetic'+x] for x in cart]
    vecList = [elecVec, magVec]
    return vecList
cmd.extend("select_vectors",select_vectors)

#Want a worker function that automates this more
#e.g. calls loadCSV, then runs a loop of select_vectors and elec_mag
#for a particular list of indices.
