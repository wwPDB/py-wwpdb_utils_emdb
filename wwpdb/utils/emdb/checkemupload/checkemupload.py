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
import inspect


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.float32):
            return float(obj)
        elif isinstance(obj, bytes):  # Add this check
            return obj.decode()  # Convert bytes to string
        return json.JSONEncoder.default(self, obj)


class EMMap:
    """
    Class representing an Electron Microscopy (EM) map.
    Initialize an EMMap instance.
    :param path2file: Path to the map file. (str)
    Load map file and extract relevant information.
    :return: None
    :raises FileNotFoundError: If the specified file is not found.
    :raises Exception: If an error occurs while loading the file or calculating the MD5 checksum.
    **Example**::
        >>> em_map = EMMap('path/to/map.mrc')
    """

    def __init__(self, path2file):
        """
        Initialize an EMMap instance.
        This function initializes an instance of EMMap with the provided map file path. It sets several attributes to None and defines some utility functions related to floating-point comparisons using a small epsilon value. Finally, it loads the map from the specified path.
        :param path2file: Path to the map file as a string.
        :returns: None
        :rtype: None
        :raises: This function does not raise any exceptions.
        **Example**::
            >>> em_map = EMMap('path/to/mapfile.em')
            expected output
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
        self.errors = []

    def load(self):
        """
        Load map file and extract relevant information.
        This function reads a map file, extracts essential information such as header details, checksum, box size, pixel size, and coordinate reference system. The function processes the file using 'mrcfile' library and then performs data manipulation to derive the necessary parameters.
        :param self: A reference to the class instance.
        :returns: None; updates attributes of the class instance with extracted information.
        :raises FileNotFoundError: If the specified file cannot be found.
        :raises Exception: If an unexpected error occurs during the file loading process.
        **Example**::
            >>> load()  # Call the load function
            expected output
        """
        try:
            with mrcfile.open(self.file, mode="r", permissive=False) as mrc:
                self.header = mrc.header
                self.hash = self.md5_checksum(mrc.data)
                self.box_size = self.header.cella.tolist()
                self.pixel_size = mrc.voxel_size.tolist()
                self.nxyz = np.array(
                    (self.header.nx, self.header.ny, self.header.nz)
                ).tolist()
                nstarts = np.array(
                    (self.header.nxstart, self.header.nystart, self.header.nzstart)
                )
                origin = nstarts * np.array(self.pixel_size)
                end = origin + np.array(self.box_size)
                (self.nstarts, self.origin, self.end) = (
                    nstarts.tolist(),
                    origin.tolist(),
                    end.tolist(),
                )
                crs = (self.header.mapc, self.header.mapr, self.header.maps)
                crsindices = (crs.index(1), crs.index(2), crs.index(3))
                if crs != (1, 2, 3):
                    self.nxyz = np.array(
                        (
                            self.nxyz[crsindices[0]],
                            self.nxyz[crsindices[1]],
                            self.nxyz[crsindices[2]],
                        )
                    ).tolist()
                    self.nstarts = np.array(
                        (
                            self.nstarts[crsindices[0]],
                            self.nstarts[crsindices[1]],
                            self.nstarts[crsindices[2]],
                        )
                    ).tolist()
        except FileNotFoundError:
            message = f"File not found: {os.path.basename(self.file)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
        except Exception as e:
            message = f"An error occurred while loading the file: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()

    def md5_checksum(self, data=None):
        """
        Calculate the MD5 checksum of the map file.
        This function calculates the MD5 checksum of the input 'data' and returns the hex digest of the resulting MD5 hash.
        :param data: The data for which the MD5 checksum is to be calculated. It should be in bytes format.
        :returns: The MD5 checksum as a hexadecimal string.
        :rtype: str
        :raises Exception: This function may raise an exception if an error occurs while calculating the MD5 checksum.
        **Example**::
            >>> md5_checksum(b'Hello, world!')
            '65a8e27d8879283831b664bd8b7f0ad4'
        """
        try:
            hash_md5 = hashlib.md5()
            # Convert the data array to bytes and update the MD5 hash
            data_bytes = data.tobytes()
            hash_md5.update(data_bytes)
            return hash_md5.hexdigest()
        except Exception as e:
            message = f"An error occurred while calculating the MD5 checksum: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()

    def same_box_size(self, another_map):
        """
        Check if the map box has the same box size as another map.
        This function checks whether the box size of the current map matches the box size of another map by comparing their box dimensions.
        :param another_map: Another EMMap instance holding the box size to compare.
        :type another_map: EMMap
        :returns: True if the box sizes are the same, False otherwise.
        :rtype: bool
        :raises: This function may raise exceptions in case of error during box size comparison.
        **Example**::
            >>> map1 = EMMap(box_size=[5, 5, 5])
            >>> map2 = EMMap(box_size=[5, 5, 5])
            >>> map3 = EMMap(box_size=[4, 4, 4])
            >>> map1.same_box_size(map2)
            True
            >>> map1.same_box_size(map3)
            False
        """
        try:
            return all(
                abs(np.array(self.box_size) - np.array(another_map.box_size)) <= self.epsilon
            )
        except Exception as e:
            message = f"An error occurred while comparing box sizes: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
            return False

    def same_or_smaller_box_size(self, another_map):
        """
        Check if the map box has a smaller or same box size than another map.
        This function compares the box sizes of two EMMap instances and returns True if the box size of the current map is smaller or the same as the box size of the other map, and False otherwise.
        :param another_map: Another EMMap instance for comparison.
        :type another_map: EMMap
        :returns: True if the current map has a smaller or same box size than the other map, False otherwise.
        :rtype: bool
        :raises: This function can raise generic exceptions if there is an error while comparing the box sizes of the maps.
        **Example**::
            >>> map1 = EMMap()
            >>> map2 = EMMap()
            >>> map1.same_or_smaller_box_size(map2)
            True
        """
        try:
            return all(np.greater_equal(np.array(self.box_size) + self.epsilon, np.array(another_map.box_size)))
        except Exception as e:
            message = f"An error occurred while comparing box sizes: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
            return False

    def overlaps(self, another_map):
        """
        Check if the map overlaps another map.
        This function compares the boundaries of two EMMap instances to determine if they overlap within a given epsilon tolerance.
        :param another_map: Another EMMap instance to compare with.
        :type another_map: EMMap
        :returns: True if the maps overlap, False otherwise.
        :rtype: bool
        :raises: This function may raise exceptions related to comparisons of map boundaries.
        **Example**::
            >>> map1 = EMMap((0, 0), (5, 5))
            >>> map2 = EMMap((3, 3), (7, 7))
            >>> map1.overlaps(map2)
            True
        """
        try:
            (origin1, end1) = (np.array(self.origin), np.array(self.end))
            (origin2, end2) = (np.array(another_map.origin), np.array(another_map.end))
            return all(abs(origin1 - origin2) <= self.epsilon) and all(
                abs(end1 - end2) <= self.epsilon
            )
        except Exception as e:
            message = f"An error occurred while comparing map boundaries: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
            return False

    def fits_inside(self, another_map):
        """
        Check if the map is completely inside another map.
        This function compares the boundaries of two EMMap instances to determine if one map is entirely contained within the other.
        :param another_map: Another EMMap instance to compare boundaries with.
        :type another_map: EMMap
        :returns: True if the map is completely inside another_map, False otherwise.
        :rtype: bool
        :raises: This function may raise exceptions in case of errors during the boundary comparison process.
        **Example**::
            >>> map1 = EMMap(origin=[1, 1], end=[4, 4])
            >>> map2 = EMMap(origin=[0, 0], end=[5, 5])
            >>> map1.fits_inside(map2)
            True
        """
        try:
            (origin1, end1) = (np.array(self.origin), np.array(self.end))
            (origin2, end2) = (np.array(another_map.origin), np.array(another_map.end))
            return all(
                np.greater_equal(origin1, origin2 - self.epsilon)
            ) and all(
                np.less_equal(end1, end2 + self.epsilon)
            )
        except Exception as e:
            message = f"An error occurred while comparing map boundaries: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
            return False

    def same_pixel_size(self, another_map):
        """
        Check if the map has the same pixel size as another map.
        This function compares the pixel sizes of two EMMap instances to determine if they are the same within a given epsilon tolerance.
        :param another_map: Another EMMap instance to compare pixel sizes with.
        :type another_map: EMMap
        :returns: True if the pixel sizes are the same, False otherwise.
        :rtype: bool
        :raises: This function may raise exceptions if there is an error during the pixel size comparison process.
        **Example**::
            >>> map1 = EMMap(pixel_size=[1.0, 1.0])
            >>> map2 = EMMap(pixel_size=[0.5, 0.5])
            >>> map1.same_pixel_size(map2)
            False
        """
        try:
            diff = abs(np.array(self.pixel_size) - np.array(another_map.pixel_size))
            return all(diff <= self.epsilon)
        except Exception as e:
            message = f"An error occurred while comparing pixel sizes: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
            return False

    def pixel_size_is_multiple(self, another_map):
        """
        Check if the pixel size of the map is a multiple of another map.
        This function compares the pixel sizes of two EMMap instances to determine if the pixel size of the calling map is a multiple of the passed map.
        :param self: The calling EMMap instance whose pixel size is being checked.
        :param another_map: Another EMMap instance to compare pixel sizes with.
        :returns: True if the pixel size is a multiple, False otherwise.
        :rtype: bool
        :raises Exception: This function catches and handles exceptions that might occur during the pixel size comparison process.
        **Example**::
            >>> em_map_1 = EMMap(pixel_size=[1.0, 1.0])
            >>> em_map_2 = EMMap(pixel_size=[0.5, 0.5])
            >>> em_map_1.pixel_size_is_multiple(em_map_2)
            True
        """
        try:
            modulo = np.array(self.pixel_size) % np.array(another_map.pixel_size)
            return all(modulo <= self.epsilon)
        except Exception as e:
            message = f"An error occurred while comparing pixel sizes: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
            return False


class Model:
    """
    Class representing a structural model in MMCIF format.
    """
    def __init__(self, path2model):
        """
        Initialize a Model instance.
        Initialize a Model instance with the path to the MMCIF file. This function sets up the necessary attributes and retrieves the coordinates for the model structure.
        :param path2model: Path to the MMCIF file.
        :type path2model: str
        :returns: None
        :rtype: None
        :raises: In case the path provided is incorrect or the file cannot be accessed, a FileNotFoundError is raised.
        **Example**::
            >>> model_instance = Model('/path/to/model.cif')
            expected output
        """
        self.file = path2model
        self.structure = None
        self.errors = []
    
    def load(self):
        """
        Load the structural model from the MMCIF file and store the coordinates.
        This function loads the structural model from the MMCIF file and stores the coordinates in the 'structure' attribute of the Model instance. Any errors encountered during the process are logged and handled accordingly.
        :param self: An instance of the Model class.
        :returns: None
        :rtype: None
        :raises: This function may raise exceptions if the file is not found or if an error occurs during the loading process.
        **Example**::
            >>> model = Model('/path/to/model.cif')
            >>> model.load()
            >>> print(model.structure)
            [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)]
        """
        self.structure = self.get_coordinates()

    def get_coordinates(self):
        """
        Get atom coordinates from the parsed MMCIF file.
        Retrieve atom coordinates from the Multi-column Coded Information File (MMCIF) data after parsing. This function extracts x, y, z coordinates from the parsed MMCIF dictionary and creates a list of coordinate tuples. Any errors encountered during parsing or coordinate extraction are logged and handled accordingly.
        :param self: Represents an instance of the class containing the MMCIF parsing functionality.
        :returns: A list of coordinate tuples (x, y, z) extracted from the MMCIF file.
        :rtype: list
        :raises ValueError: Raised if there is an issue with parsing coordinate values.
        :raises Exception: Raised in case of unexpected errors during processing.
        **Example**::
            >>> model = Model('/path/to/model.cif')
            >>> model.load()
            >>> print(model.structure)
            [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)]
        """
        # Zhe edit 24042024
        #
        # While the coordinate file has been validated to contain floats,
        # BioPython does not deal with multiple data blocks properly.
        # The second and subsequent data blocks are misparsed.
        # Therefore, the try/except handles the misfeature - rejecting
        # "bad" data
        (mmcif_dict, error, message) = self.parse_mmcif()
        if error:
            return None
        if not mmcif_dict:
            message = "No coordinates found in the MMCIF file."
            self.errors.append(message)
            print(message)
            return None
        if not (
            "_atom_site.Cartn_x" in mmcif_dict
            and "_atom_site.Cartn_y" in mmcif_dict
            and ("_atom_site.Cartn_z" in mmcif_dict)
        ):
            message = "No coordinates found in the MMCIF file."
            self.errors.append(message)
            print(message)
            return None
        structure = []
        x = mmcif_dict["_atom_site.Cartn_x"]
        y = mmcif_dict["_atom_site.Cartn_y"]
        z = mmcif_dict["_atom_site.Cartn_z"]
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
        Parse the MMCIF file and return a dictionary.
        This function reads the MMCIF file and parses its contents into a dictionary using the BioPython library. It returns the parsed dictionary, along with any errors encountered during the parsing process.
        :param self: An instance of the Model class.
        :returns: A dictionary containing the parsed MMCIF data.
        :rtype: dict
        :raises: This function may raise exceptions if the file is not found or if an error occurs during the parsing process.
        **Example**::
            >>> parse_mmcif(self)
            ({'_atom_site.Cartn_x': [...], '_atom_site.Cartn_y': [...], '_atom_site.Cartn_z': [...]}, False, "")
        """
        mmcif_dict, error, message = None, False, ""
        try:
            mmcif_dict = MMCIF2Dict(self.file)
        except FileNotFoundError:
            error = True
            message = f"File not found: {os.path.basename(self.file)}"
            self.errors.append(message)
            print(message)
        except Exception as e:
            error = True
            message = f"An error occurred while parsing the MMCIF file: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
        return (mmcif_dict, error, message)


class Validator:
    """
    Class for validating and comparing EM maps and models.
    Initializes a Validator instance with the provided EM map, half maps, and model.
    :param em_map: Primary EM map (type: EMMap). Defaults to None if not provided.
    :param half_maps: List of half maps (type: list of EMMap). Defaults to an empty list if not provided.
    :param model: Structural model (type: Model). Defaults to None if not provided.
    :returns: Results of the checks as a dictionary.
    :rtype: dict
    :raises Exception: Raised if an error occurs during the checks.
    **Example**::
        >>> validator = Validator(em_map=em_map_instance, half_maps=[half_map1, half_map2], model=model_instance)
        >>> validator.check()
        {'half_maps_to_each_other': {...}, 'primary_map_to_half_maps': [{...}, {...}], 'map_to_model': {...}}
    """

    def __init__(self, em_map=None, half_maps=None, model=None):
        """
        Initialize a Validator instance.
        Initialize a Validator instance with the provided EM map, half maps, and a structural model if available. It sets up various attributes for further validation processes.
        :param em_map: Primary EM map of type EMMap.
        :param half_maps: List of half maps of type EMMap; defaults to an empty list.
        :param model: Structural model of type Model; defaults to None.
        :returns: None.
        :rtype: None.
        :raises: No exceptions raised.
        **Example**::
            >>> validator = Validator(em_map=em_map_instance, half_maps=[first_half_map, second_half_map], model=model_instance)
        """
        self.em_map = em_map
        self.half_maps = half_maps
        self.model = model
        self.nxyz = None
        self.nstarts = None
        self.errors = []

    def check(self):
        """
        Perform a series of checks to validate and compare EM maps and structural models.
        For maps, it compares the half maps if provided, and compares each half map to the primary map. For each pair of maps, it checks if they are identical, if they have the same box size, if one has a smaller box size than the other, if they overlap, if one fits inside the other, if they have the same pixel size, and if one has a pixel size that is a multiple of the other. For the model, it checks the number and fraction of atoms outside the primary EM map.
        :param None: This function takes no parameters.
        :returns: Results of the checks in a dictionary format.
        :rtype: dict
        :raises Exception: It raises an exception if an error occurs during the checks.
        **Example**::
            >>> check()
            {'primary_map_to_half_maps': [comparison_result1, comparison_result2],
             'map_to_model': {'num_atoms_outside': num_outside, 'fraction_atoms_outside': fraction_outside}}
        """
        result = {}
        try:
            if self.half_maps:
                if len(self.half_maps) != 2:
                    result["error"] = "Two half maps must be provided."
                elif not all(
                    (os.path.isfile(half_map.file) for half_map in self.half_maps)
                ):
                    result["error"] = "One or more half maps not found."
                else:
                    result.update(
                        {
                            "half_maps_to_each_other": self._compare_maps(
                                self.half_maps[0], self.half_maps[1]
                            )
                        }
                    )
            if self.em_map and os.path.isfile(self.em_map.file):
                if self.half_maps and all(
                    (os.path.isfile(half_map.file) for half_map in self.half_maps)
                ):
                    result.update(
                        {
                            "primary_map_to_half_maps": [
                                self._compare_maps(self.em_map, half_map)
                                for half_map in self.half_maps
                            ]
                        }
                    )
                if self.model and self.model.structure:
                    (
                        num_atoms_outside,
                        fraction_atoms_outside,
                    ) = self._get_atoms_outside()
                    result["map_to_model"] = {
                        "num_atoms_outside": num_atoms_outside,
                        "fraction_atoms_outside": fraction_atoms_outside,
                    }
        except Exception as e:
            result["error"] = "An error occurred while performing checks."
            message = f"An error occurred while performing checks: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
        return result

    def _map_matrix(self, apixs, angs):
        """
        Calculate the matrix to transform Cartesian coordinates to fractional coordinates.
        This function takes arrays of apix lengths and angles in alpha, beta, gamma order as inputs to compute a matrix representing the transformation. The matrix is determined based on the provided lengths and angles according to a specific formula.
        :param apixs: Array of lengths for x, y, and z axes.
        :param angs: Array of angles in alpha, beta, gamma order.
        :returns: Array representing the transformation matrix.
        :rtype: np.ndarray
        :raises: This function may raise general exceptions if an error occurs during the matrix calculation process.
        **Example**::
            >>> _map_matrix([10, 20, 30], [90, 90, 90])
            array([[ 0.1       , -0.        ,  0.        ],
                   [ 0.        ,  0.05      ,  0.        ],
                   [ 0.        ,  0.        ,  0.03333333]])
        """
        try:
            ang = (
                angs[0] * math.pi / 180,
                angs[1] * math.pi / 180,
                angs[2] * math.pi / 180,
            )
            insidesqrt = (
                1
                + 2 * math.cos(ang[0]) * math.cos(ang[1]) * math.cos(ang[2])
                - math.cos(ang[0]) ** 2
                - math.cos(ang[1]) ** 2
                - math.cos(ang[2]) ** 2
            )
            cellvolume = apixs[0] * apixs[1] * apixs[2] * math.sqrt(insidesqrt)
            m11 = 1 / apixs[0]
            m12 = -math.cos(ang[2]) / (apixs[0] * math.sin(ang[2]))
            m13 = (
                apixs[1]
                * apixs[2]
                * (math.cos(ang[0]) * math.cos(ang[2]) - math.cos(ang[1]))
                / (cellvolume * math.sin(ang[2]))
            )
            m21 = 0
            m22 = 1 / (apixs[1] * math.sin(ang[2]))
            m23 = (
                apixs[0]
                * apixs[2]
                * (math.cos(ang[1]) * math.cos(ang[2]) - math.cos(ang[0]))
                / (cellvolume * math.sin(ang[2]))
            )
            m31 = 0
            m32 = 0
            m33 = apixs[0] * apixs[1] * math.sin(ang[2]) / cellvolume
            prematrix = [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]
            matrix = np.asarray(prematrix)
            return matrix
        except Exception as e:
            message = f"An error occurred while calculating the matrix: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
            return None

    def _matrix_indices(self, apixs, onecoor):
        """
        Calculate the matrix-derived indices for given coordinates in relation to a map.
        This function computes the matrix-derived indices for a set of coordinates based on the provided map scale and angles. It adjusts the coordinates by taking into account the starting indices of the map.
        :param apixs: The map scale value for conversion.
        :param onecoor: The coordinates for which to calculate the indices.
        :returns: A tuple of adjusted indices relative to the map's starting point.
        :rtype: Tuple[int, int, int]
        :raises Exception: If an error occurs during the index calculation process.
        **Example**::
            >>> _matrix_indices(1.0, [10, 20, 30])
            (-5, -15, -25)
        """
        try:
            angs = [
                self.em_map.header.cellb.alpha,
                self.em_map.header.cellb.beta,
                self.em_map.header.cellb.gamma,
            ]
            matrix = self._map_matrix(apixs, angs)
            result = matrix.dot(np.asarray(onecoor))
            return (
                result[0] - self.em_map.nstarts[0],
                result[1] - self.em_map.nstarts[1],
                result[2] - self.em_map.nstarts[2],
            )
        except Exception as e:
            message = f"An error occurred while calculating the indices: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
            return None

    def _header_check(self):
        """
        Check the header information of the EM map and make necessary adjustments if needed.
        This function validates the coordinate reference system of the EM map header and modifies the map dimensions and starting points accordingly if the coordinate system is not as expected.
        :param self: The instance of the class.
        :returns: None
        :raises Exception: If an error occurs during the header check process.
        **Example**::
            >>> _header_check(instance_of_class)
        """
        try:
            crs = (
                self.em_map.header.mapc,
                self.em_map.header.mapr,
                self.em_map.header.maps,
            )
            crsindices = (crs.index(1), crs.index(2), crs.index(3))
            self.nstarts = (
                self.em_map.header.nxstart,
                self.em_map.header.nystart,
                self.em_map.header.nzstart,
            )
            self.nxyz = (
                self.em_map.header.nx,
                self.em_map.header.ny,
                self.em_map.header.nz,
            )
            if crs != (1, 2, 3):
                self.nxyz = (
                    self.nxyz[crsindices[0]],
                    self.nxyz[crsindices[1]],
                    self.nxyz[crsindices[2]],
                )
                self.nstarts = (
                    self.nstarts[crsindices[0]],
                    self.nstarts[crsindices[1]],
                    self.nstarts[crsindices[2]],
                )
        except Exception as e:
            message = f"An error occurred while checking the header: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()

    def _get_indices(self, onecoor):
        """
        Find the indices corresponding to a single atom's cubic or plane coordinates.
        This function takes in a list of atom coordinates in the (x, y, z) order and determines the corresponding 8 (cubic) or 4 (plane) indices. The calculated indices are then returned as part of a tuple containing two lists: the first list holds the cubic indices, and the second list contains the float index of the input atom.
        :param onecoor: List containing the atom coordinates in (x, y, z) order.
        :return: Tuple containing two lists of indices: the first list with 8 or 4 indices in the cubic form, and the second list with the float index of the input atom.
        :raises: Exception if an error occurs while getting the indices.
        **Example**::
            >>> _get_indices([10, 20, 30])
            ([xindex, yindex, zindex], [nxstart, nystart, nzstart])
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
            (nxstart, nystart, nzstart) = self.nstarts
            if (
                self.em_map.header.cellb.alpha
                == self.em_map.header.cellb.beta
                == self.em_map.header.cellb.gamma
                == 90.0
            ):
                zindex = (
                    float(onecoor[2] - self.em_map.header.origin.z) / z_apix - nzstart
                )
                yindex = (
                    float(onecoor[1] - self.em_map.header.origin.y) / y_apix - nystart
                )
                xindex = (
                    float(onecoor[0] - self.em_map.header.origin.x) / x_apix - nxstart
                )
            else:
                apixs = [x_apix, y_apix, z_apix]
                (xindex, yindex, zindex) = self._matrix_indices(apixs, onecoor)
            oneindex = [xindex, yindex, zindex]
            return (oneindex, self.nxyz)
        except Exception as e:
            message = f"An error occurred while getting the indices: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
            return (None, None)

    def _compare_maps(self, map1, map2):
        """
        Compare two maps and return information about their similarity.
        This function compares two map objects and provides detailed information about their key properties such as voxel size, pixel size, overlap, and box size.
        :param map1: The first map object to compare. Expected type: Map.
        :param map2: The second map object to compare. Expected type: Map.
        :returns: A dictionary containing information about the comparison result, including 'em_volumes' with names and checksums of the maps, 'map_checks' with comparisons like identical hash, box size, pixel size, overlap, and box fitting.
        :rtype: dict
        :raises Exception: If an error occurs during the comparison process.
        **Example**::
            >>> map1 = Map('file1.map', 'hash1', [10, 10, 10], [0.5, 0.5, 0.5], [0, 0, 0], [9, 9, 9])
            >>> map2 = Map('file2.map', 'hash2', [10, 10, 10], [0.5, 0.5, 0.5], [0, 0, 0], [9, 9, 9])
            >>> _compare_maps(self, map1, map2)
            {'em_volumes': [{'name': 'file1.map', 'checksum': 'hash1'}, {'name': 'file2.map', 'checksum': 'hash2'}],
             'map_checks': {'identical': True, 'same_box_size': True, 'same_pixel_size': True, 'overlap': True}}
        """
        try:
            result = {
                "em_volumes": [
                    {k: v for k, v in map1.__dict__.items() if k not in ('header', 'epsilon', 'error') and not inspect.isfunction(v)},
                    {k: v for k, v in map2.__dict__.items() if k not in ('header', 'epsilon', 'error') and not inspect.isfunction(v)},
                ],
                "map_checks": {"identical": map1.hash == map2.hash},
            }
            if not result["map_checks"]["identical"]:
                same_box_size = map1.same_box_size(map2)
                result["map_checks"]["same_box_size"] = same_box_size
                if not same_box_size:
                    result["em_volumes"][0].update(
                        {"box_size": [round(x, 2) for x in map1.box_size]}
                    )
                    result["em_volumes"][1].update(
                        {"box_size": [round(x, 2) for x in map2.box_size]}
                    )
                    result["map_checks"][
                        "same_or_smaller_box_size"
                    ] = map1.same_or_smaller_box_size(map2)
                same_pixel_size = map1.same_pixel_size(map2)
                result["map_checks"]["same_pixel_size"] = same_pixel_size
                if not same_pixel_size:
                    result["em_volumes"][0].update(
                        {"pixel_size": [round(x, 2) for x in map1.pixel_size]}
                    )
                    result["em_volumes"][1].update(
                        {"pixel_size": [round(x, 2) for x in map2.pixel_size]}
                    )
                    result["map_checks"][
                        "pixel_size_is_multiple"
                    ] = map1.pixel_size_is_multiple(map2)
                overlap = map1.overlaps(map2)
                result["map_checks"]["overlap"] = overlap
                if not overlap:
                    result["em_volumes"][0].update(
                        {
                            "origin": [round(x, 2) for x in map1.origin],
                            "end": [round(x, 2) for x in map1.end],
                        }
                    )
                    result["em_volumes"][1].update(
                        {
                            "origin": [round(x, 2) for x in map2.origin],
                            "end": [round(x, 2) for x in map2.end],
                        }
                    )
                    result["map_checks"]["fits_inside"] = map1.fits_inside(map2)
            return result
        except Exception as e:
            message = f"An error occurred while comparing maps: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
            return None

    def _get_atoms_outside(self):
        """
        Get the number of atoms outside the map along with the corresponding fraction.
        This function iterates over atoms in the structure to determine the number of atoms that are outside the defined map boundaries. It calculates both the count of atoms outside the map and the fraction relative to the total number of atoms in the structure.
        :param self: The instance of the parent class containing the structure and necessary methods.
        :returns: A tuple containing the number of atoms outside the map and the fraction of such atoms relative to the total atoms.
        :rtype: Tuple[int, float]
        :raises Exception: Raised when an error occurs during the process of getting atoms outside the map.
        **Example**::
            >>> example = ClassName()
            >>> example._get_atoms_outside()
            (3, 0.05)
        """
        try:
            num_atoms_outside = 0
            for atom in self.model.structure:
                (atom_index, nxyz) = self._get_indices(atom)
                if any(
                    (
                        coord < 0 or coord > nxyz[i] - 1
                        for (i, coord) in enumerate(atom_index)
                    )
                ):
                    num_atoms_outside += 1
            fraction_atoms_outside = num_atoms_outside / len(self.model.structure)
            return (num_atoms_outside, fraction_atoms_outside)
        except Exception as e:
            message = f"An error occurred while getting atoms outside the map: {str(e)}"
            self.errors.append(message)
            print(message)
            traceback.print_exc()
            return (None, None)


def main():
    """
    Main function that parses command line arguments, performs validation checks,
    and outputs results in JSON format.
    :param None: This function does not take any parameters.
    :returns: A status code indicating success (0) or failure (1).
    :rtype: int
    :raises FileNotFoundError: If the specified input file is not found.
    :raises json.JSONDecodeError: If the specified input file is not a valid JSON file.
    :raises Exception: If an unexpected error occurs during the execution.
    **Example**::
        >>> main()
        Result written to checks.json
        0
    """
    parser = argparse.ArgumentParser(description="Perform checks on uploaded maps.")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument(
        "--input", help="Input JSON file containing paths for input files"
    )
    file_group = parser.add_argument_group(
        "file paths (required if --input is not provided)"
    )
    file_group.add_argument(
        "--primmap", help="Input MRC primary map file", required=False
    )
    file_group.add_argument(
        "--halfmaps", nargs=2, help="Input MRC half maps", required=False, default=[]
    )
    file_group.add_argument("--model", help="Input MMCIF model file", required=False)
    parser.add_argument("--output", help="Output JSON file", required=False)
    args = parser.parse_args()
    (result, errors) = ({}, [])
    # If input is provided, load the paths from the JSON file
    if args.input:
        try:
            with open(args.input, "r") as f:
                input_paths = json.load(f)
            args.primmap = input_paths.get("primmap", args.primmap)
            args.halfmaps = input_paths.get("halfmaps", args.halfmaps)
            args.model = input_paths.get("model", args.model)
        except FileNotFoundError:
            errors.append(f"Input file {args.input} not found.")
        except json.JSONDecodeError:
            errors.append(f"Input file {args.input} is not a valid JSON file.")
    if args.primmap:
        try:
            (em_map, model, half_maps) = (None, None, [])
            # Checking if files exist and loading data
            # if os.path.isfile(args.primmap):
            em_map = EMMap(args.primmap)
            em_map.load()
            errors.extend(em_map.errors)
            # Loading model if provided
            if args.model:
                model = Model(args.model)
                model.load()
                errors.extend(model.errors)
            # Loading half maps if provided
            if args.halfmaps:
                for _i, hm in enumerate(args.halfmaps, 1):
                    halfmap = EMMap(hm)
                    halfmap.load()
                    half_maps.append(halfmap)
                    errors.extend(halfmap.errors)
            validator = Validator(em_map, half_maps, model)
            result.update(validator.check())
            errors.extend(validator.errors)
        except Exception as e:
            message = f"An error occurred while performing checks: {str(e)}"
            result["error"] = message
            print(message)
            traceback.print_exc()
    elif not args.input:
        errors.append("--primmap is required when --input is not provided")
    # Adding errors to the result
    if errors:
        result["error"] = "\n".join(errors)
    # Printing results to stdout
    print(json.dumps(result, cls=NumpyEncoder))
    # Writing results to output JSON file
    if not args.output:
        if args.primmap:
            parentdir = os.path.dirname(args.primmap)
            basename = os.path.basename(args.primmap)
            (root, ext) = os.path.splitext(basename)
            while ext:
                basename = root
                (root, ext) = os.path.splitext(root)
            args.output = f"{os.path.join(parentdir, root)}-checks.json"
        else:
            args.output = "checks.json"
    with open(args.output, "w") as f:
        json.dump(result, f, indent=4, cls=NumpyEncoder)
    print(f"Result written to {args.output}")
    return 0 if "error" not in result else 1  # pylint: disable=return-in-finally,lost-exception


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as _e:  # noqa: F841
        # Use Traceback to be able to see the error in the logs
        print("An error occurred:")
        traceback.print_exc()
        sys.exit(1)
