'''
http://pymolwiki.org/index.php/cgo_arrow

(c) 2013 Thomas Holder, Schrodinger Inc.

License: BSD-2-Clause
'''

#Insert vector pair
#Default to start vector at origin of coords system
#Scaling default should be one, allow for option to scaling, increase by factor of 5
#scale separately
#If selecting atoms, shift vectors according because vectors are based around origin
#read vectors from file
#Include sphere representing radius of S tilde - optional
#Add text for actual length and angle between vectors - not scaled
#Add arc between vectors for angle

from pymol import cmd, CmdException
import pandas as pd
import argparse as ap
from pymol.cgo import *
import numpy as np
import sys

def scale_endpoint(end):
    for i in range(len(end)):
        end[i] = end[i]*4
    return end

def cgo_arrow(origin=[0,0,0], endpoint=[-3.786958749, -0.649784273, 1.59145952], radius=0.25, gap=0.0, hlength=-1, hradius=-1,
              color='blue', type='', name =''):
    '''
DESCRIPTION

    Create a CGO arrow between two picked atoms.

ARGUMENTS

    atom1 = string: single atom selection or list of 3 floats {default: sele}

    endpoint = string: list of 3 floats, currently arbitrary, normally would be desired vector from datasheet

    radius = float: arrow radius {default: 0.25}

    gap = float: gap between arrow tips and the two atoms {default: 0.0}

    hlength = float: length of head

    hradius = float: radius of head

    color = string: one or two color names {default: blue}

    name = string: name of CGO object
    '''
    from chempy import cpv
    #converting parameters to floats
    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    if type == 'electric':
        color = 'blue'
    if type == 'magnetic':
        color = 'red'
    try:
        #if they are strings, splits color1 and color2 into list of strings
        color1, color2 = color.split()
    except:
        color1 = color2 = color
    color1 = list(cmd.get_color_tuple(color1))
    color2 = list(cmd.get_color_tuple(color2))

    #get coordinates of atom1 which is 'sele'
    xyz1 = origin
    xyz2 = scale_endpoint(endpoint)
    print(xyz2)
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    #if head length parameters are not specified
    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6
    #if a gap is specified
    if gap:
        #scales normal vector by a factor of the gap
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)
    xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

    obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color1 + \
          [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color1 + color2 + \
          [1.0, 0.0]

    if not name:
        name = cmd.get_unused_name('arrow')

    cmd.load_cgo(obj, name)




data = pd.read_csv("output.csv")
print(data)
print("Select a vector to be visualized")

#cmd.extend('cgo_arrow', cgo_arrow)
