#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import math
import json
import mrcfile
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


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
        self.header = None
        self.nxyz = None
        self.nstarts = None
        self.size = None
        self.pixel_size = None
        self.epsilon = 1e-10
        self.load()

    def load(self):
        """
        Load map file and extract relevant information.

        Raises:
            FileNotFoundError: If the file is not found.
            Exception: If an error occurs during file loading.
        """
        try:
            with mrcfile.open(self.file, mode='r', permissive=True) as mrc:
                self.header = mrc.header
                self.size = [round(x, 2) for x in self.header.cella.tolist()]
                self.pixel_size = [round(x, 2) for x in mrc.voxel_size.tolist()]
                self.nxyz = np.array((self.header.nx, self.header.ny, self.header.nz)).tolist()
                self.nstarts = np.array((self.header.nxstart, self.header.nystart, self.header.nzstart)).tolist()
                crs = (self.header.mapc, self.header.mapr, self.header.maps)
                crsindices = (crs.index(1), crs.index(2), crs.index(3))
                if crs != (1, 2, 3):
                    self.nxyz = np.array((self.nxyz[crsindices[0]], self.nxyz[crsindices[1]], self.nxyz[crsindices[2]])).tolist()
                    self.nstarts = np.array((self.nstarts[crsindices[0]], self.nstarts[crsindices[1]], self.nstarts[crsindices[2]])).tolist()
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {self.file}")
        except Exception as e:
            raise Exception(f"An error occurred while loading the file: {str(e)}")

    def extremities(self):
        """
        Calculate the origin and end points of the map.

        Returns:
            tuple: Origin and end points as lists.
        """
        origin = np.array(self.nstarts) * np.array(self.pixel_size)
        end = origin + np.array(self.size)
        return origin.tolist(), end.tolist()

    def smaller_or_equal(self, another_map):
        """
        Check if the map is smaller or equal to another map in size.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if smaller or equal, False otherwise.
        """
        return all(np.array(self.size) - np.array(another_map.size) <= self.epsilon)

    def is_inside(self, another_map):
        """
        Check if the map is completely inside another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            bool: True if completely inside, False otherwise.
        """
        origin1, end1 = np.array(self.extremities())
        origin2, end2 = np.array(another_map.extremities())
        return all(origin2 - origin1 <= self.epsilon) and all(end1 - end2 <= self.epsilon)

    def acceptable_pixel_size(self, another_map):
        """
        Check if the pixel size of the map is acceptable compared to another map.

        Args:
            another_map (EMMap): Another EMMap instance.

        Returns:
            tuple: Two boolean values indicating pixel size acceptability and if it's a multiple.
        """
        diff = np.array(self.pixel_size) - np.array(another_map.pixel_size)
        modulo = np.array(self.pixel_size) % np.array(another_map.pixel_size)
        return all(diff <= self.epsilon), all(modulo <= self.epsilon)


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
        self.structure = [(float(x), float(y), float(z)) for x, y, z in
                          zip(mmcif_dict['_atom_site.Cartn_x'], mmcif_dict['_atom_site.Cartn_y'],
                              mmcif_dict['_atom_site.Cartn_z'])]


class Validator:
    """
    Class for validating and comparing EM maps and models.
    """
    def __init__(self, em_map, half_maps=[], model=None):
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

    def _map_matrix(self, apixs, angs):
        # Convert angles to radians
        angs = [math.radians(ang) for ang in angs]
        cos_alpha, cos_beta, cos_gamma = [math.cos(ang) for ang in angs]
        sin_alpha, sin_beta, sin_gamma = [math.sin(ang) for ang in angs]
        matrix = np.array([
            [apixs[0], apixs[1] * cos_gamma, apixs[2] * cos_beta],
            [0, apixs[1] * sin_gamma, -apixs[2] * sin_beta * cos_alpha],
            [0, 0, apixs[2] * sin_beta * sin_alpha]
        ])
        return matrix

    def _matrix_indices(self, apixs, onecoor):
        angs = [self.em_map.header.cellb.alpha, self.em_map.header.cellb.beta, self.em_map.header.cellb.gamma]
        matrix = self._map_matrix(apixs, angs)
        result = matrix.dot(np.asarray(onecoor))
        return result[0] - self.em_map.nstarts[0], result[1] - self.em_map.nstarts[1], result[2] - self.em_map.nstarts[2]

    def _get_indices(self, onecoor):
        apixs = [self.em_map.header.cella.x / self.em_map.nxyz[0],
                 self.em_map.header.cella.y / self.em_map.nxyz[1],
                 self.em_map.header.cella.z / self.em_map.nxyz[2]]
        oneindex = self._matrix_indices(apixs, onecoor)
        return oneindex, self.em_map.nxyz

    def check(self):
        """
        Perform a series of checks to validate and compare EM maps and structural models.
        
        For each half map, it checks whether the size is smaller or equal to the primary map,
        whether the half map is completely inside the primary map, and whether the pixel sizes are the same or multiples.

        For the model, it checks the number and fraction of atoms outside the primary EM map.
        
        Returns:
            dict: Results of the checks.
        """
        result = {
            'em_volume': os.path.basename(self.em_map.file),
            'map_checks': {},
            'model_checks': {}
        }
        
        # Check the conditions between em_map and each half_map
        for i, half_map in enumerate(self.half_maps):
            # Check if the size of half_map is smaller or equal to the primary map
            smaller_or_equal = self.em_map.smaller_or_equal(half_map)
            # Check if half_map is completely inside the primary map
            is_inside = self.em_map.is_inside(half_map)
            # Check if the pixel sizes are the same or multiples
            same_pixel_size, pixel_size_is_multiple = self.em_map.acceptable_pixel_size(half_map)             
            result['map_checks'].update({
                f'half_map_{i+1}_checks': {
                    'em_volume': os.path.basename(half_map.file),
                    'smaller_or_equal': smaller_or_equal,
                    'is_inside': is_inside,
                    'same_pixel_size': same_pixel_size,
                    'pixel_size_is_multiple': pixel_size_is_multiple
                }
            })
        
        # Check conditions for the model only if a model is provided
        if self.model:
            atoms_outside_num = 0  # Initialize counter for atoms outside the primary map
            for atom in self.model.structure:
                atom_index, nxyz = self._get_indices(atom)
                # Check if the atom is outside the primary map boundaries
                if any(coord < 0 or coord >= nxyz[i] for i, coord in enumerate(atom_index)):
                    atoms_outside_num += 1
            # Calculate the fraction of atoms outside the primary map
            atom_outside_fraction = atoms_outside_num / len(self.model.structure)
            result['model_checks'].update({
                'num_atoms_outside': atoms_outside_num,
                'fraction_atoms_outside': atom_outside_fraction
            })

        return result


def main():
    """
    Main function that parses command line arguments, performs validation checks,
    and outputs results in JSON format.
    """
    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Perform checks on uploaded maps.")
    parser.add_argument("primmap", help="Input MRC primary map file")
    parser.add_argument("--halfmaps", nargs=2, help="Input MRC half maps", required=False)
    parser.add_argument("--model", help="Input MMCIF model file", required=False)
    parser.add_argument("--output", help="Output JSON file", required=False)   
    args = parser.parse_args()

    # Checking if files exist and loading data
    if os.path.isfile(args.primmap):
        em_map = EMMap(args.primmap)
        half_maps = []
        if args.halfmaps:
            for halfmap in args.halfmaps:
                if os.path.isfile(halfmap):
                    half_maps.append(EMMap(halfmap))
        model = None
        if args.model and os.path.isfile(args.model):
            model = Model(args.model)

        # Performing validation checks
        validator = Validator(em_map, half_maps, model)
        result = validator.check()

        # Writing results to output JSON file
        parentdir = os.path.dirname(args.primmap)
        basename = os.path.basename(args.primmap)
        # Split the basename and handle double extensions
        root, ext = os.path.splitext(basename)
        if ext in ['.gz', '.bz2', '.xz']:
            root, _ = os.path.splitext(root)
        filename = args.output or f'{os.path.join(parentdir, root)}-emmapchecks.json'
        with open(filename, 'w') as f:
            json.dump(result, f, indent=4)
        print(f"Result written to {filename}")

        return 0
    return 1

# Main script execution
if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)