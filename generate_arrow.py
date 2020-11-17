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
               type='electric', name='', scaling = 7):

    from chempy import cpv
    #converting parameters to floats
    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    if type == 'electric':
        color = 'red'
        name = 'electric'+name
    if type == 'magnetic':
        color = 'blue'
        name = 'magnetic'+name
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
        length=np.linalg.norm(np.array(endpoint)-np.array(xyz1))
        xyz2 = scale_endpoint(endpoint)
        xyz2 = shift_vectors(xyz2, xyz1)
    else:
        xyz1 = origin
        length=np.linalg.norm(np.array(endpoint)-np.array(xyz1))
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


    ##Place pseudoatom with label at midpoint of vector
    v0=np.array(xyz1)
    v1=np.array(xyz2)
    loc=(v0+v1)/2

    if not name:
        name = cmd.get_unused_name('arrow')

    cmd.load_cgo(obj, f"vec_{name}")
    ##For some reason pos fails, but it just happens to put it where I want
    cmd.pseudoatom(f"lab_{name}",name="lab_"+name,label=f"{length:.2f}")#,pos=loc
    cmd.group(name,members=f"lab_{name} vec_{name}")
cmd.extend('cgo_arrow', cgo_arrow)
