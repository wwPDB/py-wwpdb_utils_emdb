#!/usr/bin/env python

import json
import mrcfile
from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import Select
import os
import sys
import argparse
import math
import numpy as np


class NotDisordered(Select):
    """

        Class used to select non-disordered atom from biopython structure instance

    """

    def accept_atom(self, atom):
        """
            Accept only atoms that at "A"
        :param atom: atom instance from biopython library
        :return: True or False
        """
        if (not atom.is_disordered()) or atom.get_altloc() == "A":
            atom.set_altloc(" ")
            return True
        else:
            return False


def hasdisorder_atom(structure):

    ress = structure.get_residues()
    disorder_flag = False
    for res in ress:
        if res.is_disordered() == 1:
            disorder_flag = True
            return disorder_flag
    return disorder_flag


def get_mapheader(map_fullname):
    header = mrcfile.open(map_fullname, permissive=True, header_only=True)
    return header.header


def get_model(model_fullname):
    """

    :param model_fullname: String of the model name
    :return: Biopython model object
    """

    p = MMCIFParser()
    io = MMCIFIO()
    orgfilename = model_fullname
    structure = p.get_structure('t', model_fullname)
    if len(structure.get_list()) > 1:
        orgmodel = model_fullname + '_org.cif'
        os.rename(model_fullname, orgmodel)
        fstructure = structure[0]
        io.set_structure(fstructure)
        io.save(model_fullname)
        usedframe = p.get_structure('first', model_fullname)
        print('!!!There are multiple models in the cif file. Here we only use the first for calculation.')
    else:
        usedframe = structure

    if hasdisorder_atom(usedframe):
        model_fullname = model_fullname + '_Alt_A.cif'
        io.set_structure(usedframe)
        print('There are alternative atom in the model here we only use A for calculations and saved as {}'
              .format(model_fullname))
        io.save(model_fullname, select=NotDisordered())
        newstructure = p.get_structure('t', model_fullname)
    else:
        newstructure = usedframe
    setattr(newstructure, "filename", orgfilename)
    tmodel = newstructure

    return tmodel


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

    # Method 2: by using the fractional coordinate matrix
    # Chosen as the main function for the current implementation

    # Figure out the order of the x, y, z based on crs info in the header
    # crs = list(self.map.header[16:19])
    # crs = [header.mapc, header.mapr, header.maps]
    # ordinds save the indices correspoding to x, y ,z
    # ordinds = [crs.index(1), crs.index(2), crs.index(3)]
    # angs = self.map.header[13:16]
    angs = [header.cellb.alpha, header.cellb.beta, header.cellb.gamma]
    matrix = map_matrix(apixs, angs)
    result = matrix.dot(np.asarray(onecoor))
    # xindex = result[0] - self.map.header[4 + ordinds[0]]
    # yindex = result[1] - self.map.header[4 + ordinds[1]]
    # zindex = result[2] - self.map.header[4 + ordinds[2]]
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
        math.cos(ang[0])**2 - \
        math.cos(ang[1])**2 - \
        math.cos(ang[2])**2

    cellvolume = apixs[0] * apixs[1] * apixs[2] * math.sqrt(insidesqrt)

    m11 = 1 / apixs[0]
    m12 = -math.cos(ang[2]) / (apixs[0] * math.sin(ang[2]))

    m13 = apixs[1] * apixs[2] * (math.cos(ang[0]) * math.cos(ang[2]) - math.cos(ang[1])) / (cellvolume * math.sin(ang[2]))
    m21 = 0
    m22 = 1 / (apixs[1] * math.sin(ang[2]))
    m23 = apixs[0] * apixs[2] * (math.cos(ang[1]) * math.cos(ang[2]) - math.cos(ang[0])) / (cellvolume * math.sin(ang[2]))
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

    # For non-cubic or skewed density maps, they might have different apix on different axises
    zdim = header.cella.z
    znintervals = header.mz
    z_apix = zdim / znintervals

    ydim = header.cella.y
    ynintervals = header.my
    y_apix = ydim / ynintervals

    xdim = header.cella.x
    xnintervals = header.mx
    x_apix = xdim / xnintervals

    # map_zsize = header.nz
    # map_ysize = header.ny
    # map_xsize = header.nx

    nxyzstart = header_check(header)
    # map_xsize, map_ysize, map_zsize = nxyzstart[0]
    nxstart, nystart, nzstart = nxyzstart[1]

###########################
    # header_check(header)
    # exit()

    # if map.header[13] == map.header[14] == map.header[15] == 90.:
    if header.cellb.alpha == header.cellb.beta == header.cellb.gamma == 90.:
        # Figure out the order of the x, y, z based on crs info in the header
        # crs = list(map.header[16:19])
        # crs = [header.mapc, header.mapr, header.maps]
        # ordinds save the indices correspoding to x, y ,z
        # ordinds = [crs.index(1), crs.index(2), crs.index(3)]

        zindex = float(onecoor[2] - header.origin.z) / z_apix - nzstart
        yindex = float(onecoor[1] - header.origin.y) / y_apix - nystart
        xindex = float(onecoor[0] - header.origin.x) / x_apix - nxstart

        # zfloor = int(floor(zindex))
        # if zfloor >= map_zsize - 1:
        #     zceil = zfloor
        # else:
        #     zceil = zfloor + 1
        #
        # yfloor = int(floor(yindex))
        # if yfloor >= map_ysize - 1:
        #     yceil = yfloor
        # else:
        #     yceil = yfloor + 1
        #
        # xfloor = int(floor(xindex))
        # if xfloor >= map_xsize - 1:
        #     xceil = xfloor
        # else:
        #     xceil = xfloor + 1
    else:
        # Method 2: by using the fractional coordinate matrix
        # Chosen as the primary for the current implementation
        apixs = [x_apix, y_apix, z_apix]
        # Method 1: by using the atom projection on planes
        # xindex, yindex, zindex = self.projection_indices(onecoor))
        xindex, yindex, zindex = matrix_indices(nxyzstart[1], apixs, onecoor, header)

    # zfloor = int(math.floor(zindex))
    # if zfloor >= map_zsize - 1:
    #     zceil = zfloor
    # else:
    #     zceil = zfloor + 1

    # yfloor = int(math.floor(yindex))
    # if yfloor >= map_ysize - 1:
    #     yceil = yfloor
    # else:
    #     yceil = yfloor + 1

    # xfloor = int(math.floor(xindex))
    # if xfloor >= map_xsize - 1:
    #     xceil = xfloor
    # else:
    #     xceil = xfloor + 1

    # indices = np.array(np.meshgrid(np.arange(xfloor, xceil + 1), np.arange(yfloor, yceil + 1),
    #                                np.arange(zfloor, zceil + 1))).T.reshape(-1, 3)
    oneindex = [xindex, yindex, zindex]

    # return (indices, oneindex)
    return oneindex, nxyzstart[0]


def check(model, header):
    atoms_outside_num = 0
    atom_counter = 0
    for atom in model.get_atoms():
        # if atom.name.startswith('H') or atom.get_parent().resname == 'HOH':
        #     continue
        atom_counter += 1
        onecoor = atom.coord
        oneindex, nxyz = get_indices(header, onecoor)
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
                     'Fraction of atoms outside box: {}'.format(os.path.basename(input_map), os.path.basename(input_model), result[0], result[1])
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
