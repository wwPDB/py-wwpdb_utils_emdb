#!/usr/bin/env python
import mrcfile
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import os
import sys
import argparse
import math
import numpy as np
import json


def get_mapheader(map_fullname):
    header = mrcfile.open(map_fullname, permissive=True, header_only=True)
    return header.header


def get_model(model_fullname):
    """

    :param model_fullname: String of the model name
    :return: Biopython model object
    """

    mmcif_dict = MMCIF2Dict(model_fullname)
    structure = []
    x = mmcif_dict['_atom_site.Cartn_x']
    y = mmcif_dict['_atom_site.Cartn_y']
    z = mmcif_dict['_atom_site.Cartn_z']
    for x, y, z in zip(x, y, z):
        try:
            newx = float(x)
            newy = float(y)
            newz = float(z)
            structure.append((newx, newy, newz))
        except ValueError:
            continue

    return structure


def header_check(header):
    crs = (header.mapc, header.mapr, header.maps)
    crsindices = (crs.index(1), crs.index(2), crs.index(3))
    nstarts = (header.nxstart, header.nystart, header.nzstart)
    nxyz = (header.nx, header.ny, header.nz)

    if crs != (1, 2, 3):
        nxyz = (nxyz[crsindices[0]], nxyz[crsindices[1]], nxyz[crsindices[2]])
        nstarts = (nstarts[crsindices[0]], nstarts[crsindices[1]], nstarts[crsindices[2]])

    return nxyz, nstarts


def matrix_indices(nstarts, apixs, onecoor, header):
    """

        Method 2: using the fractional coordinate matrix to calculate the indices when the maps are non-orthogonal

    :return:
    """

    angs = [header.cellb.alpha, header.cellb.beta, header.cellb.gamma]
    matrix = map_matrix(apixs, angs)
    result = matrix.dot(np.asarray(onecoor))
    xindex = result[0] - nstarts[0]
    yindex = result[1] - nstarts[1]
    zindex = result[2] - nstarts[2]

    return (xindex, yindex, zindex)


def map_matrix(apixs, angs):
    """

        calculate the matrix to transform Cartesian coordinates to fractional coordinates
        (check the definination to see the matrix formular)

    :param apixs: array of apix lenght
    :param angs: array of anglex in alpha, beta, gamma order
    :return:
    """

    ang = (angs[0] * math.pi / 180, angs[1] * math.pi / 180, angs[2] * math.pi / 180)
    insidesqrt = 1 + 2 * math.cos(ang[0]) * math.cos(ang[1]) * math.cos(ang[2]) - \
                 math.cos(ang[0]) ** 2 - \
                 math.cos(ang[1]) ** 2 - \
                 math.cos(ang[2]) ** 2

    cellvolume = apixs[0] * apixs[1] * apixs[2] * math.sqrt(insidesqrt)

    m11 = 1 / apixs[0]
    m12 = - math.cos(ang[2]) / (apixs[0] * math.sin(ang[2]))

    m13 = apixs[1] * apixs[2] * (math.cos(ang[0]) * math.cos(ang[2]) - math.cos(ang[1])) / (
            cellvolume * math.sin(ang[2]))
    m21 = 0
    m22 = 1 / (apixs[1] * math.sin(ang[2]))
    m23 = apixs[0] * apixs[2] * (math.cos(ang[1]) * math.cos(ang[2]) - math.cos(ang[0])) / (
            cellvolume * math.sin(ang[2]))
    m31 = 0
    m32 = 0
    m33 = apixs[0] * apixs[1] * math.sin(ang[2]) / cellvolume
    prematrix = [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]
    matrix = np.asarray(prematrix)

    return matrix


def get_indices(header, onecoor):
    """
    Find one atom's indices correspoding to its cubic or plane
    the 8 (cubic) or 4 (plane) indices are saved in indices variable
    :param map: Density map instance from TEMPy.MapParser
    :param onecoor: List contains the atom coordinates in (x, y, z) order
    :return: Tuple contains two list of index: first has the 8 or 4 indices in the cubic;
             second has the float index of the input atom
    """

    zdim = header.cella.z
    znintervals = header.mz
    z_apix = zdim / znintervals

    ydim = header.cella.y
    ynintervals = header.my
    y_apix = ydim / ynintervals

    xdim = header.cella.x
    xnintervals = header.mx
    x_apix = xdim / xnintervals

    nxyzstart = header_check(header)
    nxstart, nystart, nzstart = nxyzstart[1]

    if header.cellb.alpha == header.cellb.beta == header.cellb.gamma == 90.:
        zindex = float(onecoor[2] - header.origin.z) / z_apix - nzstart
        yindex = float(onecoor[1] - header.origin.y) / y_apix - nystart
        xindex = float(onecoor[0] - header.origin.x) / x_apix - nxstart
    else:
        apixs = [x_apix, y_apix, z_apix]
        xindex, yindex, zindex = matrix_indices(nxyzstart[1], apixs, onecoor, header)

    oneindex = [xindex, yindex, zindex]

    return oneindex, nxyzstart[0]


def check(model, header):
    atoms_outside_num = 0
    atom_counter = 0
    for atom in model:
        atom_counter += 1
        oneindex, nxyz = get_indices(header, atom)
        if oneindex[0] > nxyz[0] - 1 or oneindex[0] < 0 or \
                oneindex[1] > nxyz[1] - 1 or oneindex[1] < 0 or \
                oneindex[2] > nxyz[2] - 1 or oneindex[2] < 0:
            atoms_outside_num += 1

    atom_outside_fraction = atoms_outside_num / atom_counter if atom_counter != 0 else 0
    return atoms_outside_num, atom_outside_fraction


def read_args():
    """
        Read arguments
    :return: map, model
    """

    assert len(sys.argv) > 1, ('There has to be arguments for the commands.\n \
                              Usage main.py -m/--map  <fullmapname> -d/--model <fullmodelname>')
    parser = argparse.ArgumentParser(description='Input map and model')
    parser.add_argument('--map', '-m', nargs='?', help='EM volume')
    parser.add_argument('--model', '-d', nargs='?', help='Model')
    parser.add_argument('--output', '-o', nargs='?', help='Output JSON file')
    args = parser.parse_args()
    return args.map, args.model, args.output


def main():
    input_map, input_model, output_file = read_args()
    if os.path.isfile(input_map) and os.path.isfile(input_model):
        header = get_mapheader(input_map)
        model = get_model(input_model)
        result = check(model, header)
        result_dict = {
            'em_volume': os.path.basename(input_map),
            'model': os.path.basename(input_model),
            'atomsoutside': result[0],
            'fraction': result[1]
        }

        if output_file:
            with open(output_file, "w") as fw:
                json.dump(result_dict, fw)

        result_str = 'map:                           {}\n' \
                     'model:                         {}\n' \
                     'Number of atoms outside box:   {}\n' \
                     'Fraction of atoms outside box: {}'.format(os.path.basename(input_map),
                                                                os.path.basename(input_model), result[0], result[1])
        print(result_str)
        print('---------------------------------------')
    else:
        print('Please provide both volume map and model.')

        if output_file:
            result_dict = {
                'em_volume': '',
                'model': '',
                'atomsoutside': 0,
                'fraction': 0
            }
            with open(output_file, "w") as fw:
                json.dump(result_dict, fw)


if __name__ == '__main__':
    main()
