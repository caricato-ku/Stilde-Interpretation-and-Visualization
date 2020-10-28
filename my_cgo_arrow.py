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
from pymol.cgo import *
import numpy as np


def cgo_arrow(atom1='sele', endpoint=[2.759453589, 2.263740837, 2.505137712], radius=0.25, gap=0.0, hlength=-1, hradius=-1,
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
    if atom1=="sele":
        xyz1 = cmd.get_coords(atom1)
        xyz1 = xyz1.flatten()
        xyz1 = xyz1.tolist()
    else:
        xyz1=[float(a) for a in atom1]

    xyz2 = [float(e) for e in endpoint]
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
