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
        self.errors = []
        self.warnings = []
        self.logs = []
        self.load()

    def load(self):
        """
        Load map file and extract relevant information.
        """
        self.hash = self.md5_checksum()
        if self.hash is None:
            return

        try:
            with mrcfile.open(self.file, mode='r', permissive=False) as mrc:
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
            message = f"File not found: {os.path.basename(self.file)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
        except Exception as e:
            message = f"An error occurred while loading the file: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)

    def md5_checksum(self):
        """
        Calculate the MD5 checksum of the map file.
        """
        try:
            hash_md5 = hashlib.md5()
            with open(self.file, "rb") as f:
                # Depending on your file format, you might need to skip the header
                # before starting to read the data
                # f.seek(header_size)
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_md5.update(chunk)
            # self.hash = hash_md5.hexdigest()
            return hash_md5.hexdigest()
        except FileNotFoundError:
            message = f"File not found: {os.path.basename(self.file)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
        except Exception as e:
            message = f"An error occurred while calculating the MD5 checksum: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)

    def same_box_size(self, another_map):
        """
        Check if the map box has the same box size as another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if same box size, False otherwise.
        """
        try:
            return all(np.array(self.box_size) - np.array(another_map.box_size) <= self.epsilon)
        except Exception as e:
            message = f"An error occurred while comparing box sizes: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
            return False

    def same_or_smaller_box_size(self, another_map):
        """
        Check if the map box has a smaller or same box size than another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if smaller or same box size, False otherwise.
        """
        try:
            return all(self.le(np.array(self.box_size), np.array(another_map.box_size)))
        except Exception as e:
            message = f"An error occurred while comparing box sizes: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
            return False

    def overlaps(self, another_map):
        """
        Check if the map overlaps another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if overlaps, False otherwise.
        """
        try:
            origin1, end1 = np.array(self.origin), np.array(self.end)
            origin2, end2 = np.array(another_map.origin), np.array(another_map.end)
            # Check whether the absolute difference between the origin and end coordinates is less than epsilon
            return all(abs(origin1 - origin2) <= self.epsilon) and all(abs(end1 - end2) <= self.epsilon)
        except Exception as e:
            message = f"An error occurred while comparing map boundaries: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
            return False

    def fits_inside(self, another_map):
        """
        Check if the map is completely inside another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if completely inside, False otherwise.
        """
        try:
            origin1, end1 = np.array(self.origin), np.array(self.end)
            origin2, end2 = np.array(another_map.origin), np.array(another_map.end)
            # Check whether the origin and end coordinates are greater than or equal to the origin and end coordinates of
            # the other map, respectively
            return all(self.ge(origin1, origin2)) and all(self.le(end1, end2))
        except Exception as e:
            message = f"An error occurred while comparing map boundaries: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
            return False

    def same_pixel_size(self, another_map):
        """
        Check if the pixel size of the map is the same as another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if same pixel size, False otherwise.
        """
        try:
            diff = np.array(self.pixel_size) - np.array(another_map.pixel_size)
            return all(diff <= self.epsilon)
        except Exception as e:
            message = f"An error occurred while comparing pixel sizes: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
            return False

    def pixel_size_is_multiple(self, another_map):
        """
        Check if the pixel size of the map is a multiple of another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if pixel size is a multiple, False otherwise.
        """
        try:
            modulo = np.array(self.pixel_size) % np.array(another_map.pixel_size)
            return all(modulo <= self.epsilon)
        except Exception as e:
            message = f"An error occurred while comparing pixel sizes: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
            return False


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
        self.file = path2model
        self.errors = []
        self.warnings = []
        self.logs = []
        self.structure = self.get_coordinates()
    
    def get_coordinates(self):
        """
        Get atom coordinates from the parsed MMCIF file.
        """
        # Zhe edit 24042024
        #
        # While the coordinate file has been validated to contain floats,
        # BioPython does not deal with multiple data blocks properly.
        # The second and subsequent data blocks are misparsed.
        # Therefore, the try/except handles the misfeature - rejecting
        # "bad" data
        mmcif_dict, error, message = self.parse_mmcif()
        if error:
            return None
        if not mmcif_dict:
            message = "No coordinates found in the MMCIF file."
            self.errors.append(message)
            print(message)
            return None
        if not (
            '_atom_site.Cartn_x' in mmcif_dict
            and '_atom_site.Cartn_y' in mmcif_dict
            and '_atom_site.Cartn_z' in mmcif_dict
            ):
            message = "No coordinates found in the MMCIF file."
            self.errors.append(message)
            print(message)
            return None
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
            except ValueError as e:
                message = f"Error parsing coordinates: {str(e)}"
                self.errors.append(message)
                print(message)
                traceback.print_exc()
                continue
        return structure
    
    def parse_mmcif(self):
        """
        Parse the MMCIF file
        """
        mmcif_dict, error, message = None, False, ""
        try:
            mmcif_dict = MMCIF2Dict(self.file) # TODO: Replace with more robust parser, maybe the OneDep's one.
        except FileNotFoundError:
            error = True
            message = f"File not found: {os.path.basename(self.file)}"
            self.errors.append(message)
            print(message)
        except Exception as e:
            error = True
            message = f"An error occurred while parsing the MMCIF file: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
        return mmcif_dict, error, message


class Validator:
    """
    Class for validating and comparing EM maps and models.
    """

    def __init__(self, em_map=None, half_maps=None, model=None):
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
        self.errors = []
        self.warnings = []
        self.logs = []

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

        try:
            if self.half_maps and len(self.half_maps) == 2 and all(isinstance(half_map, EMMap) and os.path.isfile(half_map.file) for half_map in self.half_maps):
                result.update({
                    'half_maps_to_each_other': self._compare_maps(self.half_maps[0], self.half_maps[1])
                })

            if isinstance(self.em_map, EMMap) and os.path.isfile(self.em_map.file):
                if all(isinstance(half_map, EMMap) and os.path.isfile(half_map.file) for half_map in self.half_maps):
                    result.update({
                        'primary_map_to_half_maps': [self._compare_maps(self.em_map, half_map) for half_map in self.half_maps]
                    })

                if self.model and isinstance(self.model, Model) and self.model.structure:
                    num_atoms_outside, fraction_atoms_outside = self._get_atoms_outside()
                    result['map_to_model'] = {
                        'num_atoms_outside': num_atoms_outside,
                        'fraction_atoms_outside': fraction_atoms_outside
                    }

        except Exception as e:
            message = "An error occurred while performing checks." # TODO: Add more details to the error message (ask Jack Turner to help with this)
            result['warning'] = message
            message += f": {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)

        return result

    def _map_matrix(self, apixs, angs):
        """
        calculate the matrix to transform Cartesian coordinates to fractional coordinates
        (check the definination to see the matrix formular)
        :param apixs: array of apix lenght
        :param angs: array of anglex in alpha, beta, gamma order
        :return: array of matrix
        """
        try:
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
        except Exception as e:
            message = f"An error occurred while calculating the matrix: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
            return None

    def _matrix_indices(self, apixs, onecoor):
        try:
            angs = [self.em_map.header.cellb.alpha, self.em_map.header.cellb.beta, self.em_map.header.cellb.gamma]
            matrix = self._map_matrix(apixs, angs)
            result = matrix.dot(np.asarray(onecoor))
            return result[0] - self.em_map.nstarts[0], result[1] - self.em_map.nstarts[1], result[2] - self.em_map.nstarts[
                2]
        except Exception as e:
            message = f"An error occurred while calculating the indices: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
            return None

    def _header_check(self):
        try:
            crs = (self.em_map.header.mapc, self.em_map.header.mapr, self.em_map.header.maps)
            crsindices = (crs.index(1), crs.index(2), crs.index(3))
            self.nstarts = (self.em_map.header.nxstart, self.em_map.header.nystart, self.em_map.header.nzstart)
            self.nxyz = (self.em_map.header.nx, self.em_map.header.ny, self.em_map.header.nz)

            if crs != (1, 2, 3):
                self.nxyz = (self.nxyz[crsindices[0]], self.nxyz[crsindices[1]], self.nxyz[crsindices[2]])
                self.nstarts = (self.nstarts[crsindices[0]], self.nstarts[crsindices[1]], self.nstarts[crsindices[2]])
        except Exception as e:
            message = f"An error occurred while checking the header: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)

    def _get_indices(self, onecoor):
        """
            Find one atom's indices correspoding to its cubic or plane
            the 8 (cubic) or 4 (plane) indices are saved in indices variable
        :param onecoor: List contains the atom coordinates in (x, y, z) order
        :return: Tuple contains two list of index: first has the 8 or 4 indices in the cubic;second has the float index of the input atom
        """
        try:
            zdim = self.em_map.header.cella.z
            znintervals = self.em_map.header.mz
            z_apix = zdim / znintervals

            ydim = self.em_map.header.cella.y
            ynintervals = self.em_map.header.my
            y_apix = ydim / ynintervals

            xdim = self.em_map.header.cella.x
            xnintervals = self.em_map.header.mx
            x_apix = xdim / xnintervals

            self._header_check()
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
        except Exception as e:
            message = f"An error occurred while getting the indices: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
            return None, None

    def _compare_maps(self, map1, map2):
        try:
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
        except Exception as e:
            message = f"An error occurred while comparing maps: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
            return None

    def _get_atoms_outside(self):
        try:
            num_atoms_outside = 0  # Initialize counter for atoms outside the primary map
            for atom in self.model.structure:
                atom_index, nxyz = self._get_indices(atom)
                # Check if the atom is outside the primary map boundaries
                if any(coord < 0 or coord > nxyz[i] - 1 for i, coord in enumerate(atom_index)):
                    num_atoms_outside += 1
            # Calculate the fraction of atoms outside the primary map
            fraction_atoms_outside = num_atoms_outside / len(self.model.structure)
            return num_atoms_outside, fraction_atoms_outside
        except Exception as e:
            message = f"An error occurred while getting atoms outside the map: {str(e)}\n"
            message += traceback.format_exc()
            self.logs.append(message)
            print(message)
            return None, None

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

    result, errors, warnings, logs = {}, [], [], []
    try:
        em_map, model, half_maps = None, None, []
        # Loading data
        em_map = EMMap(args.primmap)
        errors.extend(em_map.errors)
        warnings.extend(em_map.warnings)
        logs.extend(em_map.logs)

        # Loading model if provided
        if args.model:
            model = Model(args.model)
            errors.extend(model.errors)
            warnings.extend(model.warnings)
            logs.extend(model.logs)
        
        # Loading half maps if provided
        if args.halfmaps:
            for i, hm in enumerate(args.halfmaps, 1):
                halfmap = EMMap(hm)
                half_maps.append(halfmap)
                errors.extend(halfmap.errors)
                warnings.extend(halfmap.warnings)
                logs.extend(halfmap.logs)
        
        # Performing validation checks
        validator = Validator(em_map, half_maps, model)
        result.update(validator.check())
        errors.extend(validator.errors)
        warnings.extend(validator.warnings)
        logs.extend(validator.logs)

        # Adding errors, warnings, and logs to the result
        if errors:
            result['error'] = '\n'.join(errors)
        if warnings:
            result['warning'] = '\n'.join(warnings)
        if logs:
            result['log'] = '\n'.join(logs)

    except Exception as e:
        message = "An error occurred while performing checks"
        result['warning'] = message
        message += f": {str(e)}\n"
        message += traceback.format_exc()
        result['log'] = message
        print(message)

    finally:
        parentdir = os.path.dirname(args.primmap)
        basename = os.path.basename(args.primmap)
        root, ext = os.path.splitext(basename)
        while ext:
            basename = root
            root, ext = os.path.splitext(root)
        filename = args.output or f'{os.path.join(parentdir, root)}-checks.json'
        message = f"\nResults written to {filename}"
        result['log'] += message
        print(message)
        # Printing results to stdout
        print(json.dumps(result, indent=4))
        # Writing results to output JSON file
        with open(filename, 'w') as f:
            json.dump(result, f, indent=4)
        return 0 if 'error' not in result else 1
        
# Main script execution
if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as _e:  # noqa: F841
        # Use Traceback to be able to see the error in the logs
        print(f"An error occurred: {_e}")
        traceback.print_exc()
        sys.exit(1)
