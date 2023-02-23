#!/usr/bin/env python3

import argparse
import mrcfile
import numpy as np
import json
import sys


def le(array1, array2) -> bool:
    result = list(map(lambda x: x[0] <= x[1], list(zip(array1, array2))))
    return sum(result) == len(result)


def multiple(array1, array2):
    result = list(map(lambda x: x[0] % x[1] == 0, list(zip(array1, array2))))
    return sum(result) == len(result)


class LoadMap:

    def __init__(self, path2file):
        self.file = path2file
        self.dimensions = None
        self.size = None
        self.offset = None
        self.pixel_size = None

    def load(self):
        with mrcfile.open(self.file, 'r') as mrc:
            self.dimensions = np.array((mrc.header.nx, mrc.header.ny, mrc.header.nz)).tolist()
            self.size = mrc.header.cella.tolist()
            self.offset = np.array((mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart)).tolist()
            self.pixel_size = mrc.voxel_size.tolist()

    def extremities(self):
        origin = np.array(self.offset) * np.array(self.pixel_size)
        end = origin + np.array(self.size)
        return origin.tolist(), end.tolist()

    def smaller_or_equal(self, another_map):
        return le(self.size, another_map.size)

    def is_inside(self, another_map):
        origin1, end1 = self.extremities()
        origin2, end2 = another_map.extremities()
        return le(origin2, origin1) and le(end1, end2)

    # Needs to be split into two new functions
    def accepted_pixel_size(self, another_map):
        pxlsz1 = list(map(lambda x: round(x, 2), self.pixel_size))
        pxlsz2 = list(map(lambda x: round(x, 2), another_map.pixel_size))
        if not multiple(pxlsz1, pxlsz2):
            print(f"WARNING! Pixel sizes of first map ({pxlsz1}) are not multiple of second map`s sizes ({pxlsz2}).")
        return le(pxlsz2, pxlsz1)


class UploadMapCheck:

    def __init__(self, json_data):
        self.input = json_data
        self.output = None

    def check_all_maps(self):
        primmap = LoadMap(self.input['primary_map'])
        primmap.load()
        halfmaps = [LoadMap(self.input['half_maps'][i]) for i in range(len(self.input['half_maps']))]
        for halfmap in halfmaps:
            halfmap.load()

        self.output = {
            'primary_map': {
                'path': primmap.file,
                'dimensions': primmap.dimensions,
                'size': primmap.size,
                'offset': primmap.offset,
                'pixel_size': primmap.pixel_size,
                'smaller_or_equal': {},
                'is_inside': {},
                'accepted_pixel_size': {}
            },
            'half_maps': {}
        }
        for halfmap in halfmaps:
            self.output['primary_map']['smaller_or_equal'][halfmap.file] = primmap.smaller_or_equal(halfmap)
            self.output['primary_map']['is_inside'][halfmap.file] = primmap.is_inside(halfmap)
            self.output['primary_map']['accepted_pixel_size'][halfmap.file] = primmap.accepted_pixel_size(halfmap)
            self.output['half_maps'][halfmap.file] = {
                'dimensions': halfmap.dimensions,
                'size': halfmap.size,
                'offset': halfmap.offset,
                'pixel_size': halfmap.pixel_size
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
    Takes the paths from a json file, then do some checks on maps 
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
