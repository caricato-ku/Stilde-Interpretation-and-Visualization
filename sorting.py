    ###Add comments to describe what subsections do

"""Top level description of sorting.py here"""

if _name_=="_main_":
    import numpy as np
    import pandas as pd
    import argparse as ap
    from pathlib import Path



    ##Read file as mandatory argument, error if not given.
    ##Read output as optional with -o flag, "output.csv" is default
    ##Can get description of input by calling with -h
    parser=ap.ArgumentParser()
    parser.add_argument('file',help="CSV file of stilde values")
    parser.add_argument('-o','--output',help='Sorted stilde output file, with extension',
                        default="output.csv")
    parser.add_argument('-n','--number',help="Number of stilde rows to sort",type=int,default=20)
    args=parser.parse_args()

    filename = args.file
    output = args.output
    number = args.number


    try:
        dataframe = pd.read_csv(filename)
    except FileNotFoundError as err:
        print("Exception: {0}".format(err))
        exit()
    print('\n')


    dataframe.insert(7, 'Electric Magnitude', np.nan, False)
    dataframe.insert(10, 'Magnetic Magnitude', np.nan, False)
    dataframe.insert(11, 'Cosine of Angle', np.nan, False)


    ##Numpy functions can be quickly/succinctly applied across dataframe
    elecLength=np.linalg.norm(dataframe[['ElectricX','ElectricY','ElectricZ']].values,axis=1)
    magLength=np.linalg.norm(dataframe[['MagneticX','MagneticY','MagneticZ']].values,axis=1)
    dot=np.einsum('ij,ij->i',dataframe[['ElectricX','ElectricY','ElectricZ']].values,
               dataframe[['MagneticX','MagneticY','MagneticZ']].values)

    dataframe['Electric Magnitude']=elecLength
    dataframe['Magnetic Magnitude']=magLength
    dataframe['Cosine of Angle']=dot/(elecLength*magLength)


    ##Removes need for Negative? column
    dataframe['absS']= dataframe['S'].apply(abs)
    dataframe.sort_values(by='absS',axis=0,inplace=True,ascending = False)
    dataframe.drop('absS', axis=1, inplace=True)
    ##Requires pandas 1.1.0
    #dataframe.sort_values(by='S',axis=0,inplace=True,ascending = False,key=np.abs)
    print(dataframe)


    number = int(input("Enter the number of S-tilde values you'd like to have in the exported file\n"))

    newData = dataframe[:number]
    print(f"Top {number} s-tilde values")
    print(newData)


    ##Path will automatically enter address compatible with Linux/Mac/Windows
    outdir = Path('./SortedData')
    if not Path.exists(outdir):
        Path.mkdir(outdir)
    ##Can extend Path objects with /
    fullname = outdir / output


    try:
        ##Removed row label
        newData.to_csv(fullname,index=False)
        print("New file saved successfully")
    except OSError as err:
        print("Directory error: {0}".format(err))

    ##Should try to break the program into functions and a main section
    ##Google "main guarding python"

    #visual representation of s-tilde
    #print a relative number of s-tilde based on percent error/total sum of all values
