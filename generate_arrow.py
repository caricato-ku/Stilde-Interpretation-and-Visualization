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
from pymol.cgo import *

def scale_endpoint(end, factor=5):
    for i in range(len(end)):
        end[i] = end[i]*factor
    return end

def shift_vectors(start, atomCoords):
    for i in range(len(start)):
        start[i] += atomCoords[i]
    return start

def cgo_arrow(origin, endpoint, color='blue', radius=0.25, gap=0.0, hlength=-1,  hradius=-1,
               type='electric', name=''):

    from chempy import cpv
    #converting parameters to floats
    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    if type == 'electric':
        color = 'red'
        name = 'electric'
    if type == 'magnetic':
        color = 'blue'
        name = 'magnetic'
    try:
        color1, color2 = color.split()
    except:
        color1 = color2 = color
    color1 = list(cmd.get_color_tuple(color1))
    color2 = list(cmd.get_color_tuple(color2))

    if origin == 'sele':
        xyz1 = cmd.get_coords('sele', 1)
        xyz1 = xyz1.flatten()
        xyz1 = xyz1.tolist()
        xyz2 = scale_endpoint(endpoint)
        xyz2 = shift_vectors(xyz2, xyz1)
    else:
        xyz1 = origin
        xyz2 = scale_endpoint(endpoint)

    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6
    if gap:
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
