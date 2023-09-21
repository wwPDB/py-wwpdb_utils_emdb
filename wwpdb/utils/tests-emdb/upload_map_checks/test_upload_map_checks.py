import unittest
from wwpdb.utils.emdb.upload_map_checks.upload_map_checks import *
import numpy as np

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

mock_map1 = {
    'path': '../data/em/mock_map1.map',
    'dimensions': [100, 100, 100],
    'offset': [50, 50, 50],
    'pixel_size': [1., 1., 1.]
}
mock_map2 = {
    'path': '../data/em/mock_map2.map',
    'dimensions': [200, 200, 200],
    'offset': [0, 0, 0],
    'pixel_size': [.75, .75, .75]
}
mock_map3 = {
    'path': '../data/em/mock_map3.map',
    'dimensions': [200, 200, 200],
    'offset': [0, 0, 0],
    'pixel_size': [.5, .5, .5]
}

mock_mrc_file(**mock_map1)
mock_mrc_file(**mock_map2)
mock_mrc_file(**mock_map3)


class TestLoadMap(unittest.TestCase):

    def setUp(self):
        self.load_map1 = LoadMap(mock_map1['path'])
        self.load_map1.load()
        self.load_map2 = LoadMap(mock_map2['path'])
        self.load_map2.load()
        self.load_map3 = LoadMap(mock_map3['path'])
        self.load_map3.load()

    def test_load_sets_attributes_correctly(self):
        self.assertEqual(self.load_map1.dimensions, [100, 100, 100])
        self.assertEqual(self.load_map1.size, [100., 100., 100.])
        self.assertEqual(self.load_map1.offset, [50, 50, 50])
        self.assertEqual(self.load_map1.pixel_size, [1., 1., 1.])

    def test_extremities(self):
        origin, end = self.load_map1.extremities()
        self.assertEqual(origin, [50., 50., 50.])
        self.assertEqual(end, [150., 150., 150.])

    def test_smaller_or_equal_returns_true_when_other_map_is_smaller_or_equal(self):
        self.assertTrue(self.load_map1.smaller_or_equal(self.load_map1))
        self.assertTrue(self.load_map1.smaller_or_equal(self.load_map2))
        self.assertTrue(self.load_map1.smaller_or_equal(self.load_map3))
        self.assertTrue(self.load_map2.smaller_or_equal(self.load_map2))
        self.assertTrue(self.load_map3.smaller_or_equal(self.load_map1))
        self.assertTrue(self.load_map3.smaller_or_equal(self.load_map2))

    def test_smaller_or_equal_returns_false_when_other_map_is_bigger(self):
        self.assertFalse(self.load_map2.smaller_or_equal(self.load_map1))
        self.assertFalse(self.load_map2.smaller_or_equal(self.load_map3))

    def test_is_inside_returns_true_when_other_map_contains_current_map(self):
        self.assertTrue(self.load_map1.is_inside(self.load_map1))
        self.assertTrue(self.load_map1.is_inside(self.load_map2))
        self.assertTrue(self.load_map2.is_inside(self.load_map2))
        self.assertTrue(self.load_map3.is_inside(self.load_map2))

    def test_is_inside_returns_false_when_other_map_do_not_contains_current_map(self):
        self.assertFalse(self.load_map1.is_inside(self.load_map3))
        self.assertFalse(self.load_map2.is_inside(self.load_map1))
        self.assertFalse(self.load_map2.is_inside(self.load_map3))
        self.assertFalse(self.load_map3.is_inside(self.load_map1))
    def test_acceptable_pixel_size_returns_true_when_other_map_pixel_size_is_small_enough(self):
        self.assertTrue(self.load_map1.acceptable_pixel_size(self.load_map1)[0])
        self.assertTrue(self.load_map1.acceptable_pixel_size(self.load_map2)[0])
        self.assertTrue(self.load_map1.acceptable_pixel_size(self.load_map3)[0])
        self.assertTrue(self.load_map2.acceptable_pixel_size(self.load_map2)[0])
        self.assertTrue(self.load_map2.acceptable_pixel_size(self.load_map3)[0])

    def test_acceptable_pixel_size_returns_false_when_other_map_pixel_size_is_not_small_enough(self):
        self.assertFalse(self.load_map2.acceptable_pixel_size(self.load_map1)[0])
        self.assertFalse(self.load_map3.acceptable_pixel_size(self.load_map1)[0])
        self.assertFalse(self.load_map3.acceptable_pixel_size(self.load_map2)[0])

    def test_acceptable_pixel_size_returns_true_when_other_map_pixel_size_is_multiple(self):
        self.assertTrue(self.load_map1.acceptable_pixel_size(self.load_map1)[1])
        self.assertTrue(self.load_map1.acceptable_pixel_size(self.load_map3)[1])
        self.assertTrue(self.load_map2.acceptable_pixel_size(self.load_map2)[1])


    def test_acceptable_pixel_size_returns_false_when_other_map_pixel_size_is_not_multiple(self):
        self.assertFalse(self.load_map1.acceptable_pixel_size(self.load_map2)[1])
        self.assertFalse(self.load_map2.acceptable_pixel_size(self.load_map1)[1])
        self.assertFalse(self.load_map2.acceptable_pixel_size(self.load_map3)[1])
        self.assertFalse(self.load_map3.acceptable_pixel_size(self.load_map1)[1])
        self.assertFalse(self.load_map3.acceptable_pixel_size(self.load_map2)[1])


class TestUploadMapCheck(unittest.TestCase):
    def setUp(self):
        self.json_data1 = {
            "primary_map": mock_map1['path'],
            "other_maps": [
                mock_map2['path'],
                mock_map3['path']
            ]
        }
        self.checker1 = UploadMapCheck(self.json_data1)
        self.checker1.check_all_maps()

        self.json_data2 = {
            "primary_map": mock_map2['path'],
            "other_maps": [
                mock_map1['path'],
                mock_map3['path']
            ]
        }
        self.checker2 = UploadMapCheck(self.json_data2)
        self.checker2.check_all_maps()

        self.json_data3 = {
            "primary_map": mock_map3['path'],
            "other_maps": [
                mock_map1['path'],
                mock_map2['path']
            ]
        }
        self.checker3 = UploadMapCheck(self.json_data3)
        self.checker3.check_all_maps()

    def test_primary_map(self):
        self.assertEqual(self.json_data1['primary_map'], self.checker1.output['primary_map']['path'])
        self.assertEqual(mock_map1['dimensions'], self.checker1.output['primary_map']['dimensions'])
        self.assertEqual(
            (np.array(mock_map1['dimensions']) * np.array(mock_map1['pixel_size'])).tolist(),
            self.checker1.output['primary_map']['size']
        )
        self.assertEqual(mock_map1['offset'], self.checker1.output['primary_map']['offset'])
        self.assertEqual(mock_map1['pixel_size'], self.checker1.output['primary_map']['pixel_size'])

        self.assertEqual(self.json_data2['primary_map'], self.checker2.output['primary_map']['path'])
        self.assertEqual(mock_map2['dimensions'], self.checker2.output['primary_map']['dimensions'])
        self.assertEqual(
            (np.array(mock_map2['dimensions']) * np.array(mock_map2['pixel_size'])).tolist(),
            self.checker2.output['primary_map']['size']
        )
        self.assertEqual(mock_map2['offset'], self.checker2.output['primary_map']['offset'])
        self.assertEqual(self.checker2.output['primary_map']['pixel_size'], mock_map2['pixel_size'])

        self.assertEqual(self.json_data3['primary_map'], self.checker3.output['primary_map']['path'])
        self.assertEqual(mock_map3['dimensions'], self.checker3.output['primary_map']['dimensions'])
        self.assertEqual(
            (np.array(mock_map3['dimensions']) * np.array(mock_map3['pixel_size'])).tolist(),
            self.checker3.output['primary_map']['size']
        )
        self.assertEqual(mock_map3['offset'], self.checker3.output['primary_map']['offset'])
        self.assertEqual(mock_map3['pixel_size'], self.checker3.output['primary_map']['pixel_size'])

    def test_other_maps(self):
        checks = {
            'smaller_or_equal': [
                [True, True],
                [False, False],
                [True, True]
            ],
            'is_inside': [
                [True, False],
                [False, False],
                [False, True]
            ],
            'acceptable_pixel_size': [
                [True, True],
                [False, True],
                [False, False]
            ],
            'pixel_size_is_multiple': [
                [False, True],
                [False, False],
                [False, False]
            ]
        }
        for i, checker in enumerate([self.checker1, self.checker2, self.checker3]):
            for j, other_map in enumerate(checker.output['other_maps']):
                for check in checks:
                    self.assertEqual(checks[check][i][j], checker.output['primary_map'][check][other_map])

if __name__ == '__main__':
    unittest.main()
