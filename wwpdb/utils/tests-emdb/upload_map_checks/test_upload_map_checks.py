import unittest
from unittest.mock import patch
from wwpdb.utils.emdb.upload_map_checks.upload_map_checks import *
import numpy as np


def mock_load():
    mock_file_path = '/path/to/mock/file'

    map1, map2, map3 = LoadMap(mock_file_path), LoadMap(mock_file_path), LoadMap(mock_file_path)

    map1.dimensions = (100, 100, 100)
    map1.size = (100, 100, 100)
    map1.offset = (50, 50, 50)
    map1.pixel_size = (1., 1., 1.)

    map2.dimensions = (200, 200, 200)
    map2.size = (150, 150, 150)
    map2.offset = (0, 0, 0)
    map2.pixel_size = (.75, .75, .75)

    map3.dimensions = (200, 200, 200)
    map3.size = (100, 100, 100)
    map3.offset = (0, 0, 0)
    map3.pixel_size = (.5, .5, .5)

    return map1, map2, map3


def mock_mrc_file(
        path='/path/to/mock/file.map',
        dimensions=(10, 10, 10),
        offset=(0, 0, 0),
        pixel_size=(1., 1., 1.)
    ):
    """

    :param path:
    :param dimensions:
    :param offset:
    :param pixel_size:
    :return:
    """
    aux = 1
    for d in dimensions:
        aux *= d
    mock_data = np.zeros(aux, dtype=np.int8).reshape(dimensions)
    with mrcfile.new(path, overwrite=True) as mrc:
        mrc.set_data(mock_data)
        mrc.header['nx'], mrc.header['ny'], mrc.header['nz'] = dimensions
        mrc.header['nxstart'], mrc.header['nystart'], mrc.header['nzstart'] = offset
        mrc.voxel_size = pixel_size
        mrc.update_header_from_data()
        mrc.update_header_stats()
    return True

def set_input(primmap, halfmap1, halfmap2):
    mock_mrc_file(**primmap)
    mock_mrc_file(**halfmap1)
    mock_mrc_file(**halfmap2)
    data = {
        'primary_map': primmap['path'],
        'other_maps': [
            halfmap1['path'],
            halfmap2['path']
        ]
    }
    return data

def other_maps(halfmap1, halfmap2):
    data = {
        halfmap1['path']: {
            'dimensions': halfmap1['dimensions'],
            'size': [
                round(x, 2) for x in (
                        np.array(halfmap1['dimensions']) * np.array(halfmap1['pixel_size'])
                ).tolist()
            ],
            'offset': halfmap1['offset'],
            'pixel_size': halfmap1['pixel_size']
        },
        halfmap2['path']: {
            'dimensions': halfmap2['dimensions'],
            'size': [
                round(x, 2) for x in (
                        np.array(halfmap2['dimensions']) * np.array(halfmap2['pixel_size'])
                ).tolist()
            ],
            'offset': halfmap2['offset'],
            'pixel_size': halfmap2['pixel_size']
        }
    }
    return data

mock_map1 = {
    'path': '../data/em/mock_map1.map',
    'dimensions': (100, 100, 100),
    'offset': (50, 50, 50),
    'pixel_size': (1., 1., 1.)
}
mock_map2 = {
    'path': '../data/em/mock_map2.map',
    'dimensions': (200, 200, 200),
    'offset': (0, 0, 0),
    'pixel_size': (.75, .75, .75)
}
mock_map3 = {
    'path': '../data/em/mock_map3.map',
    'dimensions': (200, 200, 200),
    'offset': (0, 0, 0),
    'pixel_size': (.5, .5, .5)
}
mock_primmap_A = {
    'path': '../data/em/mock_primmap_A.map',
    'dimensions': [240, 240, 240],
    'offset': [0, 0, 0],
    'pixel_size': [1.08, 1.08, 1.08]
}
mock_halfmap_A1 = {
    'path': '../data/em/mock_halfmap_A1.map',
    'dimensions': [240, 240, 240],
    'offset': [0, 0, 0],
    'pixel_size': [1.08, 1.08, 1.08]
}
mock_halfmap_A2 = {
    'path': '../data/em/mock_halfmap_A2.map',
    'dimensions': [240, 240, 240],
    'offset': [0, 0, 0],
    'pixel_size': [1.08, 1.08, 1.08]
}
mock_primmap_B = {
    'path': '../data/em/mock_primmap_B.map',
    'dimensions': [350, 350, 350],
    'offset': [0, 0, 0],
    'pixel_size': [0.86, 0.86, 0.86]
}


class MyTestCase(unittest.TestCase):

    def test_le_false(self):
        self.assertEqual(le((1., 1., 1.), (0., 1., 2.)), False)

    def test_le_true(self):
        self.assertEqual(le((1., 1., 1.), (1., 1., 2.)), True)

    def test_multiple_false(self):
        self.assertEqual(multiple((1., 1., 1.), (0., 1., 2.)), False)

    def test_multiple_true(self):
        self.assertEqual(multiple((6., 6., 6.), (3., 3., 3.)), True)

    @patch.object(LoadMap, 'extremities')
    def test_extremities(self, mock_extremities):
        map1, _, _ = mock_load()
        mock_extremities.return_value = (50, 150)
        self.assertEqual(map1.extremities(), (50, 150))

    @patch.object(LoadMap, 'smaller_or_equal')
    def test_smaller_or_equal_true(self, mock_smaller_or_equal):
        map1, map2, _ = mock_load()
        mock_smaller_or_equal.return_value = True
        self.assertEqual(map1.smaller_or_equal(map2), True)

    @patch.object(LoadMap, 'smaller_or_equal')
    def test_smaller_or_equal_false(self, mock_smaller_or_equal):
        map1, map2, _ = mock_load()
        mock_smaller_or_equal.return_value = False
        self.assertEqual(map2.smaller_or_equal(map1), False)

    @patch.object(LoadMap, 'is_inside')
    def test_is_inside_true(self, mock_is_inside):
        map1, map2, _ = mock_load()
        mock_is_inside.return_value = True
        self.assertEqual(map1.is_inside(map2), True)

    @patch.object(LoadMap, 'is_inside')
    def test_is_inside_false(self, mock_is_inside):
        map1, map2, _ = mock_load()
        mock_is_inside.return_value = False
        self.assertEqual(map2.is_inside(map1), False)

    @patch.object(LoadMap, 'acceptable_pixel_size')
    def test_acceptable_pixel_size_true(self, mock_acceptable_pixel_size):
        map1, map2, _ = mock_load()
        mock_acceptable_pixel_size.return_value = (True, False)
        self.assertEqual(map1.acceptable_pixel_size(map2), (True, False))

    @patch.object(LoadMap, 'acceptable_pixel_size')
    def test_acceptable_pixel_size_false(self, mock_acceptable_pixel_size):
        map1, map2, _ = mock_load()
        mock_acceptable_pixel_size.return_value = (False, False)
        self.assertEqual(map2.acceptable_pixel_size(map1), (False, False))

    @patch.object(LoadMap, 'acceptable_pixel_size')
    def test_acceptable_pixel_size_true_multiple(self, mock_acceptable_pixel_size):
        map1, _, map2 = mock_load()
        mock_acceptable_pixel_size.return_value = (True, True)
        self.assertEqual(map1.acceptable_pixel_size(map2), (True, True))

    def test_load(self):
        mock_mrc_file(**mock_primmap_A)
        mock_map = LoadMap(mock_primmap_A['path'])
        mock_map.load()
        test = (
            mock_map.file,
            mock_map.dimensions,
            mock_map.offset,
            mock_map.pixel_size
        )
        target = (
            mock_primmap_A['path'],
            mock_primmap_A['dimensions'],
            mock_primmap_A['offset'],
            mock_primmap_A['pixel_size']
        )
        self.assertEqual(test, target)

    def test_check_all_maps_A_A1_A2(self):
        data = set_input(mock_primmap_A, mock_halfmap_A1, mock_halfmap_A2)
        mock_check = UploadMapCheck(data)
        test = mock_check.check_all_maps()
        target = {
            'primary_map': {
                'path': mock_primmap_A['path'],
                'dimensions': mock_primmap_A['dimensions'],
                'size': [
                    round(x,2) for x in (
                            np.array(mock_primmap_A['dimensions']) * np.array(mock_primmap_A['pixel_size'])
                    ).tolist()
                ],
                'offset': mock_primmap_A['offset'],
                'pixel_size': mock_primmap_A['pixel_size'],
                'smaller_or_equal': {
                    mock_halfmap_A1['path']: True,
                    mock_halfmap_A2['path']: True
                },
                'is_inside': {
                    mock_halfmap_A1['path']: True,
                    mock_halfmap_A2['path']: True
                },
                'acceptable_pixel_size': {
                    mock_halfmap_A1['path']: True,
                    mock_halfmap_A2['path']: True
                },
                'pixel_size_is_multiple': {
                    mock_halfmap_A1['path']: True,
                    mock_halfmap_A2['path']: True
                }
            },
            'other_maps': other_maps(mock_halfmap_A1, mock_halfmap_A2)
        }
        self.assertEqual(test, target)

    def test_check_all_maps_B_A1_A2(self):
        data = set_input(mock_primmap_B, mock_halfmap_A1, mock_halfmap_A2)
        mock_check = UploadMapCheck(data)
        test = mock_check.check_all_maps()
        target = {
            'primary_map': {
                'path': mock_primmap_B['path'],
                'dimensions': mock_primmap_B['dimensions'],
                'size': [
                    round(x, 2) for x in (
                            np.array(mock_primmap_B['dimensions']) * np.array(mock_primmap_B['pixel_size'])
                    ).tolist()
                ],
                'offset': mock_primmap_B['offset'],
                'pixel_size': mock_primmap_B['pixel_size'],
                'smaller_or_equal': {
                    mock_halfmap_A1['path']: False,
                    mock_halfmap_A2['path']: False
                },
                'is_inside': {
                    mock_halfmap_A1['path']: False,
                    mock_halfmap_A2['path']: False
                },
                'acceptable_pixel_size': {
                    mock_halfmap_A1['path']: False,
                    mock_halfmap_A2['path']: False
                },
                'pixel_size_is_multiple': {
                    mock_halfmap_A1['path']: False,
                    mock_halfmap_A2['path']: False
                }
            },
            'other_maps': other_maps(mock_halfmap_A1, mock_halfmap_A2)
        }
        self.assertEqual(test, target)


if __name__ == '__main__':
    unittest.main()
