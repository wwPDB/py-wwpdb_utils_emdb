#!/usr/bin/env python3

import argparse
import mrcfile
import numpy as np
import json
import sys


class LoadMap:

    def __init__(self, path2file):
        self.file = path2file
        self.dimensions = None
        self.size = None
        self.offset = None
        self.pixel_size = None

    def load(self):
        with mrcfile.open(self.file, mode='r', permissive=True) as mrc:
            self.dimensions = np.array((mrc.header.nx, mrc.header.ny, mrc.header.nz)).tolist()
            self.size = [round(x, 2) for x in mrc.header.cella.tolist()]
            self.offset = np.array((mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart)).tolist()
            self.pixel_size = [round(x, 2) for x in mrc.voxel_size.tolist()]

    def extremities(self):
        origin = np.array(self.offset) * np.array(self.pixel_size)
        end = origin + np.array(self.size)
        return origin.tolist(), end.tolist()

    def smaller_or_equal(self, another_map):
        return all(self.size <= another_map.size)

    def is_inside(self, another_map):
        origin1, end1 = self.extremities()
        origin2, end2 = another_map.extremities()
        return all(origin2 <= origin1) and all(end1 <= end2)

    def acceptable_pixel_size(self, another_map):
        return all(another_map.pixel_size, self.pixel_size), all(self.pixel_size % another_map.pixel_size == 0)


class UploadMapCheck:

    def __init__(self, json_data):
        self.input = json_data
        self.output = None

    def check_all_maps(self):
        primmap = LoadMap(self.input['primary_map'])
        primmap.load()
        othermaps = [LoadMap(self.input['other_maps'][i]) for i in range(len(self.input['other_maps']))]
        for othermap in othermaps:
            othermap.load()

        self.output = {
            'primary_map': {
                'path': primmap.file,
                'dimensions': primmap.dimensions,
                'size': primmap.size,
                'offset': primmap.offset,
                'pixel_size': primmap.pixel_size,
                'smaller_or_equal': {},
                'is_inside': {},
                'acceptable_pixel_size': {},
                'pixel_size_is_multiple':{}
            },
            'other_maps': {}
        }
        for othermap in othermaps:
            self.output['primary_map']['smaller_or_equal'][othermap.file] = primmap.smaller_or_equal(othermap)
            self.output['primary_map']['is_inside'][othermap.file] = primmap.is_inside(othermap)
            self.output['primary_map']['acceptable_pixel_size'][othermap.file], \
                self.output['primary_map']['pixel_size_is_multiple'][othermap.file] = primmap.acceptable_pixel_size(othermap)
            self.output['other_maps'][othermap.file] = {
                'dimensions': othermap.dimensions,
                'size': othermap.size,
                'offset': othermap.offset,
                'pixel_size': othermap.pixel_size
            }
        return self.output

    def __repr__(self):
        if not self.output:
            return str(self.input)
        return str(self.output)


def run_uploaded_map_checks(input_json_file, output_json_file):
    with open(input_json_file, 'r') as infile:
        data = json.load(infile)
    result = UploadMapCheck(data).check_all_maps()
    with open(output_json_file, 'w') as outfile:
        json.dump(result, outfile, separators=(',', ':'), sort_keys=False, indent=4)
    return True


def main():
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
    sys.exit(0) if main() else sys.exit(1)
