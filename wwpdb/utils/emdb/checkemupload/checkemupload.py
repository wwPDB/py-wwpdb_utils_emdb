#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import math
import json
import mrcfile
import hashlib
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import traceback


class EMMap:
    """
    Class representing an Electron Microscopy (EM) map.
    """

    def __init__(self, path2file):
        """
        Initialize an EMMap instance.

        Args:
            path2file (str): Path to the map file.
        """
        self.file = path2file
        self.hash = None
        self.header = None
        self.nxyz = None
        self.nstarts = None
        self.box_size = None
        self.pixel_size = None
        self.origin = None
        self.end = None
        # Epsilon to use as tolerance for floating point comparisons
        self.epsilon = np.finfo(np.float32).eps
        # Lambda to check whether a value is greater than or equal to another value, by using epsilon as tolerance
        self.ge = lambda x, y: x >= y - self.epsilon
        # Lambda to check whether a value is less than or equal to another value, by using epsilon as tolerance
        self.le = lambda x, y: x <= y + self.epsilon
        self.load()

    def load(self):
        """
        Load map file and extract relevant information.

        Raises:
            FileNotFoundError: If the file is not found.
            Exception: If an error occurs during file loading.
        """
        try:
            self.md5_checksum()
            with mrcfile.mmap(self.file, mode='r', permissive=False) as mrc:
                self.header = mrc.header
                self.box_size = self.header.cella.tolist()
                self.pixel_size = mrc.voxel_size.tolist()
                self.nxyz = np.array((self.header.nx, self.header.ny, self.header.nz)).tolist()
                nstarts = np.array((self.header.nxstart, self.header.nystart, self.header.nzstart))
                origin = nstarts * np.array(self.pixel_size)
                end = origin + np.array(self.box_size)
                self.nstarts, self.origin, self.end = nstarts.tolist(), origin.tolist(), end.tolist()
                crs = (self.header.mapc, self.header.mapr, self.header.maps)
                crsindices = (crs.index(1), crs.index(2), crs.index(3))
                if crs != (1, 2, 3):
                    self.nxyz = np.array((
                        self.nxyz[crsindices[0]],
                        self.nxyz[crsindices[1]],
                        self.nxyz[crsindices[2]]
                    )).tolist()
                    self.nstarts = np.array((
                        self.nstarts[crsindices[0]],
                        self.nstarts[crsindices[1]],
                        self.nstarts[crsindices[2]]
                    )).tolist()
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {self.file}")  # pylint: disable=raise-missing-from
        except Exception as e:
            raise Exception(  # pylint: disable=broad-exception-raised,raise-missing-from
                f"An error occurred while loading the file: {str(e)}")

    def md5_checksum(self):
        """
        Calculate the MD5 checksum of the map file.
        """
        hash_md5 = hashlib.md5()
        with open(self.file, "rb") as f:
            # Depending on your file format, you might need to skip the header
            # before starting to read the data
            # f.seek(header_size)
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        self.hash = hash_md5.hexdigest()

    def same_box_size(self, another_map):
        """
        Check if the map box has the same box size as another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if same box size, False otherwise.
        """
        return all(np.array(self.box_size) - np.array(another_map.box_size) <= self.epsilon)

    def same_or_smaller_box_size(self, another_map):
        """
        Check if the map box has a smaller or same box size than another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if smaller or same box size, False otherwise.
        """
        return all(self.le(np.array(self.box_size), np.array(another_map.box_size)))

    def overlaps(self, another_map):
        """
        Check if the map overlaps another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if overlaps, False otherwise.
        """
        origin1, end1 = np.array(self.origin), np.array(self.end)
        origin2, end2 = np.array(another_map.origin), np.array(another_map.end)
        # Check whether the absolute difference between the origin and end coordinates is less than epsilon
        return all(abs(origin1 - origin2) <= self.epsilon) and all(abs(end1 - end2) <= self.epsilon)

    def fits_inside(self, another_map):
        """
        Check if the map is completely inside another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if completely inside, False otherwise.
        """
        origin1, end1 = np.array(self.origin), np.array(self.end)
        origin2, end2 = np.array(another_map.origin), np.array(another_map.end)
        # Check whether the origin and end coordinates are greater than or equal to the origin and end coordinates of
        # the other map, respectively
        return all(self.ge(origin1, origin2)) and all(self.le(end1, end2))

    def same_pixel_size(self, another_map):
        """
        Check if the pixel size of the map is the same as another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if same pixel size, False otherwise.
        """
        diff = np.array(self.pixel_size) - np.array(another_map.pixel_size)
        return all(diff <= self.epsilon)

    def pixel_size_is_multiple(self, another_map):
        """
        Check if the pixel size of the map is a multiple of another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if pixel size is a multiple, False otherwise.
        """
        modulo = np.array(self.pixel_size) % np.array(another_map.pixel_size)
        return all(modulo <= self.epsilon)


class Model:
    """
    Class representing a structural model in MMCIF format.
    """

    def __init__(self, path2model):
        """
        Initialize a Model instance.

        Args:
            path2model (str): Path to the MMCIF file.
        """
        # Extracting atom coordinates from the MMCIF file
        mmcif_dict = MMCIF2Dict(path2model)
        self.file = path2model
        # Zhe edit 24042024
        #
        # While the coordinate file has been validated to contain floats,
        # BioPython does not deal with multiple data blocks properly.
        # The second and subsequent data blocks are misparsed.
        # Therefore, the try/except handles the misfeature - rejecting
        # "bad" data
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
        self.structure = structure


class Validator:
    """
    Class for validating and comparing EM maps and models.
    """

    def __init__(self, em_map, half_maps=None, model=None):
        """
        Initialize a Validator instance.

        Args:
            em_map (EMMap): Primary EM map.
            half_maps (list of EMMap, optional): List of half maps. Defaults to [].
            model (Model, optional): Structural model. Defaults to None.
        """
        self.em_map = em_map
        self.half_maps = half_maps
        self.model = model
        self.nxyz = None
        self.nstarts = None

    def _map_matrix(self, apixs, angs):
        """
        calculate the matrix to transform Cartesian coordinates to fractional coordinates
        (check the definination to see the matrix formular)
        :param apixs: array of apix lenght
        :param angs: array of anglex in alpha, beta, gamma order
        :return: array of matrix
        """

        ang = (angs[0] * math.pi / 180, angs[1] * math.pi / 180, angs[2] * math.pi / 180)
        insidesqrt = (1 + 2 * math.cos(ang[0]) * math.cos(ang[1]) * math.cos(ang[2]) - math.cos(ang[0]) ** 2 - math.cos(ang[1]) ** 2 - math.cos(ang[2]) ** 2)
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

    def _matrix_indices(self, apixs, onecoor):
        angs = [self.em_map.header.cellb.alpha, self.em_map.header.cellb.beta, self.em_map.header.cellb.gamma]
        matrix = self._map_matrix(apixs, angs)
        result = matrix.dot(np.asarray(onecoor))
        return result[0] - self.em_map.nstarts[0], result[1] - self.em_map.nstarts[1], result[2] - self.em_map.nstarts[
            2]

    def header_check(self):
        crs = (self.em_map.header.mapc, self.em_map.header.mapr, self.em_map.header.maps)
        crsindices = (crs.index(1), crs.index(2), crs.index(3))
        self.nstarts = (self.em_map.header.nxstart, self.em_map.header.nystart, self.em_map.header.nzstart)
        self.nxyz = (self.em_map.header.nx, self.em_map.header.ny, self.em_map.header.nz)

        if crs != (1, 2, 3):
            self.nxyz = (self.nxyz[crsindices[0]], self.nxyz[crsindices[1]], self.nxyz[crsindices[2]])
            self.nstarts = (self.nstarts[crsindices[0]], self.nstarts[crsindices[1]], self.nstarts[crsindices[2]])

    def _get_indices(self, onecoor):
        """
                Find one atom's indices correspoding to its cubic or plane
                the 8 (cubic) or 4 (plane) indices are saved in indices variable
            :param onecoor: List contains the atom coordinates in (x, y, z) order
            :return: Tuple contains two list of index: first has the 8 or 4 indices in the cubic;second has the float index of the input atom
            """

        zdim = self.em_map.header.cella.z
        znintervals = self.em_map.header.mz
        z_apix = zdim / znintervals

        ydim = self.em_map.header.cella.y
        ynintervals = self.em_map.header.my
        y_apix = ydim / ynintervals

        xdim = self.em_map.header.cella.x
        xnintervals = self.em_map.header.mx
        x_apix = xdim / xnintervals

        self.header_check()
        # map_xsize, map_ysize, map_zsize = self.nxyz
        nxstart, nystart, nzstart = self.nstarts

        if self.em_map.header.cellb.alpha == self.em_map.header.cellb.beta == self.em_map.header.cellb.gamma == 90.:
            # crs = [self.em_map.header.mapc, self.em_map.header.mapr, self.em_map.header.maps]
            # ordinds = [crs.index(1), crs.index(2), crs.index(3)]
            zindex = float(onecoor[2] - self.em_map.header.origin.z) / z_apix - nzstart
            yindex = float(onecoor[1] - self.em_map.header.origin.y) / y_apix - nystart
            xindex = float(onecoor[0] - self.em_map.header.origin.x) / x_apix - nxstart
        else:
            apixs = [x_apix, y_apix, z_apix]
            xindex, yindex, zindex = self._matrix_indices(apixs, onecoor)

        oneindex = [xindex, yindex, zindex]

        return oneindex, self.nxyz

    def _compare_maps(self, map1, map2):
        result = {
            'em_volumes': [
                {
                    'name': os.path.basename(map1.file),
                    'checksum': map1.hash
                },
                {
                    'name': os.path.basename(map2.file),
                    'checksum': map2.hash
                }
            ],
            'map_checks': {
                'identical': map1.hash == map2.hash
            }
        }
        if not result['map_checks']['identical']:
            # Check if the maps have the same box size
            same_box_size = map1.same_box_size(map2)
            result['map_checks']['same_box_size'] = same_box_size
            if not same_box_size:
                result['em_volumes'][0].update({
                    'box_size': [round(x, 2) for x in map1.box_size]
                })
                result['em_volumes'][1].update({
                    'box_size': [round(x, 2) for x in map2.box_size]
                })
                # Check if the map box has a smaller or same box size than another map
                result['map_checks']['same_or_smaller_box_size'] = map1.same_or_smaller_box_size(map2)
            # Check if the maps have the same pixel size
            same_pixel_size = map1.same_pixel_size(map2)
            result['map_checks']['same_pixel_size'] = same_pixel_size
            if not same_pixel_size:
                result['em_volumes'][0].update({
                    'pixel_size': [round(x, 2) for x in map1.pixel_size]
                })
                result['em_volumes'][1].update({
                    'pixel_size': [round(x, 2) for x in map2.pixel_size]
                })
                # Check if the maps have pixel size that is a multiple
                result['map_checks']['pixel_size_is_multiple'] = map1.pixel_size_is_multiple(map2)
            # Check if the maps overlap
            overlap = map1.overlaps(map2)
            result['map_checks']['overlap'] = overlap
            if not overlap:
                result['em_volumes'][0].update({
                    'origin': [round(x, 2) for x in map1.origin],
                    'end': [round(x, 2) for x in map1.end]
                })
                result['em_volumes'][1].update({
                    'origin': [round(x, 2) for x in map2.origin],
                    'end': [round(x, 2) for x in map2.end]
                })
                # Check if the maps fit inside each other
                result['map_checks']['fits_inside'] = map1.fits_inside(map2)

        return result

    def _get_atoms_outside(self):
        num_atoms_outside = 0  # Initialize counter for atoms outside the primary map
        for atom in self.model.structure:
            atom_index, nxyz = self._get_indices(atom)
            # Check if the atom is outside the primary map boundaries
            if any(coord < 0 or coord > nxyz[i] - 1 for i, coord in enumerate(atom_index)):
                num_atoms_outside += 1
        # Calculate the fraction of atoms outside the primary map
        fraction_atoms_outside = num_atoms_outside / len(self.model.structure)
        return num_atoms_outside, fraction_atoms_outside

    def check(self):
        """
        Perform a series of checks to validate and compare EM maps and structural models.

        For maps, it compares the half maps if provided, and compares each half map to the primary map.
        For each pair of maps, it checks if they are identical, if they have the same box size, if one has a smaller box
        size than the other, if they overlap, if one fits inside the other, if they have the same pixel size, and if
        one has a pixel size that is a multiple of the other.

        For the model, it checks the number and fraction of atoms outside the primary EM map.

        Returns:
            dict: Results of the checks.
        """
        result = {}

        if self.half_maps:
            if len(self.half_maps) != 2:
                raise Exception("Two half maps must be provided.")  # pylint: disable=broad-exception-raised
            if not all(os.path.isfile(half_map.file) for half_map in self.half_maps):
                raise FileNotFoundError("One or more half maps not found.")
            result.update({
                'half_maps_to_each_other': self._compare_maps(self.half_maps[0], self.half_maps[1]),
                # Compare half maps to each other
                'primary_map_to_half_maps': [self._compare_maps(self.em_map, half_map) for half_map in self.half_maps]
                # Compare primary map to each half map
            })

        if self.model and self.model.structure:
            num_atoms_outside, fraction_atoms_outside = self._get_atoms_outside()
            result['map_to_model'] = {
                'num_atoms_outside': num_atoms_outside,
                # Calculate the fraction of atoms outside the primary map
                'fraction_atoms_outside': fraction_atoms_outside
            }

        return result


def main():
    """
    Main function that parses command line arguments, performs validation checks,
    and outputs results in JSON format.
    """
    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Perform checks on uploaded maps.")
    parser.add_argument("primmap", help="Input MRC primary map file")
    parser.add_argument("--halfmaps", nargs=2, help="Input MRC half maps", required=False, default=[])
    parser.add_argument("--model", help="Input MMCIF model file", required=False)
    parser.add_argument("--output", help="Output JSON file", required=False)
    args = parser.parse_args()

    # Checking if files exist and loading data
    if os.path.isfile(args.primmap):
        em_map = EMMap(args.primmap)
        half_maps = []
        for halfmap in args.halfmaps:
            if os.path.isfile(halfmap):
                half_maps.append(EMMap(halfmap))
        model = None
        if args.model and os.path.isfile(args.model):
            model = Model(args.model)

        # Performing validation checks
        validator = Validator(em_map, half_maps, model)
        result = validator.check()

        # Printing results to stdout
        print(json.dumps(result, indent=4))

        # Writing results to output JSON file
        parentdir = os.path.dirname(args.primmap)
        basename = os.path.basename(args.primmap)
        root, ext = os.path.splitext(basename)
        while ext:
            basename = root
            root, ext = os.path.splitext(root)
        filename = args.output or f'{os.path.join(parentdir, root)}-checks.json'
        with open(filename, 'w') as f:
            json.dump(result, f, indent=4)
        print(f"Result written to {filename}")

        return 0
    return 1


# Main script execution
if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as _e:  # noqa: F841
        # Use Traceback to be able to see the error in the logs
        print("An error occurred:")
        traceback.print_exc()
        sys.exit(1)
