#!/usr/bin/env python3

import argparse
import mrcfile
import numpy as np
import json
import sys


class LoadMap:
    """
    Class for loading and processing map files.
    """

    def __init__(self, path2file):
        """
        Initialize a LoadMap instance.

        Args:
            path2file (str): Path to the map file.
        """
        self.file = path2file
        self.dimensions = None
        self.size = None
        self.offset = None
        self.pixel_size = None
        self.epsilon = 1e-10

    def load(self):
        """
        Load map file and extract relevant information.

        Raises:
            FileNotFoundError: If the file is not found.
            Exception: If an error occurs during file loading.
        """
        try:
            with mrcfile.open(self.file, mode='r', permissive=True) as mrc:
                self.dimensions = np.array((mrc.header.nx, mrc.header.ny, mrc.header.nz)).tolist()
                self.size = [round(x, 2) for x in mrc.header.cella.tolist()]
                self.offset = np.array((mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart)).tolist()
                self.pixel_size = [round(x, 2) for x in mrc.voxel_size.tolist()]
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
        origin = np.array(self.offset) * np.array(self.pixel_size)
        end = origin + np.array(self.size)
        return origin.tolist(), end.tolist()

    def smaller_or_equal(self, another_map):
        """
        Check if the map is smaller or equal to another map in size.

        Args:
            another_map (LoadMap): Another LoadMap instance.

        Returns:
            bool: True if smaller or equal, False otherwise.
        """
        return all(np.array(self.size) - np.array(another_map.size) <= self.epsilon)

    def is_inside(self, another_map):
        """
        Check if the map is completely inside another map.

        Args:
            another_map (LoadMap): Another LoadMap instance.

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
            another_map (LoadMap): Another LoadMap instance.

        Returns:
            tuple: Two boolean values indicating pixel size acceptability and if it's a multiple.
        """
        diff = np.array(self.pixel_size) - np.array(another_map.pixel_size)
        modulo = np.array(self.pixel_size) % np.array(another_map.pixel_size)
        return all(diff <= self.epsilon), all(modulo <= self.epsilon)

    def acceptable_pixel_size(self, another_map):
        """
        Check if the pixel size of the map is acceptable compared to another map.

        Args:
            another_map (LoadMap): Another LoadMap instance.

        Returns:
            tuple: Two boolean values indicating pixel size acceptability and if it's a multiple.
        """
        return all(another_map.pixel_size <= self.pixel_size), all(self.pixel_size % another_map.pixel_size == 0)


class UploadMapCheck:
    """
    Class for checking uploaded maps.
    """

    def __init__(self, json_data):
        """
        Initialize an UploadMapCheck instance.

        Args:
            json_data (dict): JSON data containing map information.
        """
        self.input_data = json_data
        self.output_data = None

    def check_all_maps(self):
        """
        Check all maps provided in the JSON input.

        Returns:
            dict: Dictionary containing the check results.
        """
        primmap = LoadMap(self.input_data.get('primary_map', None))
        primmap.load()
        othermaps = [LoadMap(self.input_data['other_maps'][i]) for i in range(len(self.input_data['other_maps']))]
        for othermap in othermaps:
            othermap.load()

        self.output_data = {
            'primary_map': {
                'path': primmap.file,
                'dimensions': primmap.dimensions,
                'size': primmap.size,
                'offset': primmap.offset,
                'pixel_size': primmap.pixel_size,
                'smaller_or_equal': {},
                'is_inside': {},
                'acceptable_pixel_size': {},
                'pixel_size_is_multiple': {}
            },
            'other_maps': {}
        }
        for othermap in othermaps:
            self.output_data['primary_map']['smaller_or_equal'][othermap.file] = primmap.smaller_or_equal(othermap)
            self.output_data['primary_map']['is_inside'][othermap.file] = primmap.is_inside(othermap)
            self.output_data['primary_map']['acceptable_pixel_size'][othermap.file], \
                self.output_data['primary_map']['pixel_size_is_multiple'][othermap.file] = primmap.acceptable_pixel_size(othermap)
            self.output_data['other_maps'][othermap.file] = {
                'dimensions': othermap.dimensions,
                'size': othermap.size,
                'offset': othermap.offset,
                'pixel_size': othermap.pixel_size
            }
        return self.output_data

    def __repr__(self):
        """
        Return a string representation of the UploadMapCheck instance.

        Returns:
            str: String representation.
        """
        if not self.output_data:
            return str(self.input_data)
        return str(self.output_data)


def run_uploaded_map_checks(input_json_file, output_json_file):
    """
    Run the map checks and write the results to an output JSON file.

    Args:
        input_json_file (str): Path to the input JSON file.
        output_json_file (str): Path to the output JSON file.

    Returns:
        bool: True if the checks were successful, False otherwise.
    """
    try:
        with open(input_json_file, 'r') as infile:
            data = json.load(infile)
        result = UploadMapCheck(data).check_all_maps()
        with open(output_json_file, 'w') as outfile:
            json.dump(result, outfile, separators=(',', ':'), sort_keys=False, indent=4)
        return True
    except FileNotFoundError:
        print(f"Error: Input JSON file not found: {input_json_file}")
        return False
    except Exception as e:
        print(f"Error: {str(e)}")
        return False


def main():
    """
    Main function to parse command-line arguments and execute map checks.

    Returns:
        int: 0 if successful, 1 if there was an error.
    """
    description = """
    Takes the path from a json file, then do some checks on maps 
    and writes the results in another json file.
    """
    usage = """
    upload_map_checks.py -i[--input] <input json file> -o[--output] <output json file>
    """
    parser = argparse.ArgumentParser(description=description, usage=usage)
    parser.add_argument("-i", "--input", type=str, dest="input", required=True,
                        help="Path to  INPUT json file.")
    parser.add_argument("-o", "--output", type=str, dest="output", required=True,
                        help="Path to  OUTPUT json file.")
    args = parser.parse_args()
    return run_uploaded_map_checks(args.input, args.output)


if __name__ == "__main__":
    try:
        result = main()
        sys.exit(0 if result else 1)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)
