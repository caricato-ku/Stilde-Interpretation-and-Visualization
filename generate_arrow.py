'''
http://pymolwiki.org/index.php/cgo_arrow

(c) 2013 Thomas Holder, Schrodinger Inc.

License: BSD-2-Clause
'''


#If selecting atoms, shift vectors according because vectors are based around origin
#read vectors from file
#Include sphere representing radius of S tilde - optional
#Add text for actual length and angle between vectors - not scaled
#Add arc between vectors for angle

from pymol import cmd, CmdException
import argparse as ap
from pymol.cgo import *
import numpy as np
import sys

def scale_endpoint(end, factor=5):
    print(end)
    for i in range(len(end)):
        end[i] = end[i]*factor
    return end

def cgo_arrow(origin=[0,0,0], endpoint=[], radius=0.25, gap=0.0, hlength=-1, hradius=-1,
              color='blue', type='', name =''):

    from chempy import cpv
    #converting parameters to floats
    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    if type == 'electric':
        color = 'red'
    if type == 'magnetic':
        color = 'blue'
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
cmd.extend('cgo_arrow', cgo_arrow)
