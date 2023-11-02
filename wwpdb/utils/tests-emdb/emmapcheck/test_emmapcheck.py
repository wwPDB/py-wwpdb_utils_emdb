import unittest
from unittest.mock import mock_open, patch, MagicMock
import numpy as np
from wwpdb.utils.emdb.emmapcheck.emmapcheck import *

# Example content of a dummy MMCIF file
DUMMY_CIF_CONTENT = """
data_test
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
...
ATOM   1   N N   . ALA A 1 1   ? 11.104 22.040 12.034 1.00 10.00 ? 1  ALA A N   1 
ATOM   2   C CA  . ALA A 1 1   ? 11.000 20.634 12.034 1.00 10.00 ? 1  ALA A CA  1 
...
"""

class TestEMMap(unittest.TestCase):

    def setUp(self):
        # Define the mock map properties
        self.mock_map_properties = [
            {
                'path': '../data/em/mock_map1.map',
                'nxyz': [100, 100, 100],
                'nstarts': [50, 50, 50],
                'pixel_size': [1., 1., 1.]
            },
            {
                'path': '../data/em/mock_map2.map',
                'nxyz': [200, 200, 200],
                'nstarts': [0, 0, 0],
                'pixel_size': [.75, .75, .75]
            },
            {
                'path': '../data/em/mock_map3.map',
                'nxyz': [200, 200, 200],
                'nstarts': [0, 0, 0],
                'pixel_size': [.5, .5, .5]
            }
        ]

        # Create mock maps based on the properties
        self.mock_maps = [self.create_mock_map(props) for props in self.mock_map_properties]

        # Load the maps into EMMap objects
        self.load_map1, self.load_map2, self.load_map3 = [self.create_em_map(mock_map) for mock_map in self.mock_maps]

    def create_mock_map(self, props):
        mock_mrc = Mock()
        mock_mrc.data.shape = props['nxyz']
        mock_mrc.header.nxstart, mock_mrc.header.nystart, mock_mrc.header.nzstart = props['nstarts']
        mock_mrc.voxel_size = Mock()
        mock_mrc.voxel_size.flags = Mock(x=props['pixel_size'][0], y=props['pixel_size'][1], z=props['pixel_size'][2])
        return mock_mrc

    @patch("mrcfile.open")
    def create_em_map(self, mock_mrc, mock_open):
        mock_open.return_value.__enter__.return_value = mock_mrc
        return EMMap(mock_mrc.header.path)

    @patch("mrcfile.open")
    def load_mock_map(self, path, mock_open):
        if path == '../data/em/mock_map1.map':
            mock_open.return_value.__enter__.return_value = self.mock_mrc1
        elif path == '../data/em/mock_map2.map':
            mock_open.return_value.__enter__.return_value = self.mock_mrc2
        elif path == '../data/em/mock_map3.map':
            mock_open.return_value.__enter__.return_value = self.mock_mrc3
        return EMMap(path)

    def test_smaller_or_equal_returns_true_when_other_map_is_smaller_or_equal(self):
        load_map1 = self.load_mock_map('../data/em/mock_map1.map')
        load_map2 = self.load_mock_map('../data/em/mock_map2.map')
        load_map3 = self.load_mock_map('../data/em/mock_map3.map') 
        
    def test_extremities(self):
        origin, end = self.em_map.extremities()
        self.assertTrue(all(isinstance(coord, float) for coord in origin))
        self.assertTrue(all(isinstance(coord, float) for coord in end))

    def test_smaller_or_equal(self):
        # Comparing against the same map, should be True
        self.assertTrue(self.em_map.smaller_or_equal(self.em_map))

    def test_is_inside(self):
        # Comparing against the same map, should be True
        self.assertTrue(self.em_map.is_inside(self.em_map))

    def test_acceptable_pixel_size(self):
        # Comparing against the same map, should be True for both checks
        same_pixel_size, pixel_size_is_multiple = self.em_map.acceptable_pixel_size(self.em_map)
        self.assertTrue(same_pixel_size)
        self.assertTrue(pixel_size_is_multiple)

    def test_load_sets_attributes_correctly(self):
        self.assertEqual(self.load_map1.nxyz, [100, 100, 100])
        self.assertEqual(self.load_map1.size, [100., 100., 100.])
        self.assertEqual(self.load_map1.nstarts, [50, 50, 50])
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


class TestModel(unittest.TestCase):

    def test_init(self):
        # Assuming a dummy MMCIF file is created with atom coordinates
        model = Model('path/to/dummy.cif')
        self.assertIsNotNone(model.structure)
        # Assuming the structure is a list of tuples (x, y, z)
        self.assertTrue(all(isinstance(coord, tuple) and len(coord) == 3 for coord in model.structure))


class TestValidator(unittest.TestCase):

    def setUp(self):
        self.em_map = EMMap(mock_map1['path'])
        self.half_maps = [EMMap(mock_map2['path']), EMMap(mock_map3['path'])]
        self.model = Model('path/to/dummy.cif')
        self.validator = Validator(self.em_map, self.half_maps, self.model)

    @patch('your_module.Validator._map_matrix', return_value=np.array([[1,0,0],[0,1,0],[0,0,1]]))
    @patch('your_module.Validator._matrix_indices', return_value=(1, 1, 1))
    def test_check(self, mock_matrix_indices, mock_map_matrix):
        result = self.validator.check()

        # Check the presence of keys
        self.assertIn('em_volume', result)
        self.assertIn('map_checks', result)
        self.assertIn('model_checks', result)

        # Further assertions can be added to check the content of the results

    def test_map_matrix(self):
        matrix = self.validator._map_matrix([1,1,1], [90,90,90])
        self.assertTrue(isinstance(matrix, np.ndarray))

    def test_matrix_indices(self):
        indices = self.validator._matrix_indices([1,1,1], [1,1,1])
        self.assertTrue(all(isinstance(coord, float) for coord in indices))

    def test_get_indices(self):
        indices, nxyz = self.validator._get_indices([1,1,1])
        self.assertTrue(all(isinstance(coord, float) for coord in indices))
        self.assertTrue(all(isinstance(coord, int) for coord in nxyz))


class TestEndToEnd(unittest.TestCase):

    def test_main(self):
        # Mock command line arguments and call the main function
        sys.argv = ['your_module.py', 'path/to/dummy.mrc', '--halfmaps', 'path/to/dummy_half1.mrc', 'path/to/dummy_half2.mrc', '--model', 'path/to/dummy.cif']
        main()

        # Check that the output file is created and contains the expected results
        with open('expected_output.json', 'r') as f:
            expected_output = json.load(f)

        self.assertTrue(os.path.exists('path/to/dummy-emmapchecks.json'))
        with open('path/to/dummy-emmapchecks.json', 'r') as f:
            actual_output = json.load(f)

        self.assertEqual(expected_output, actual_output)

if __name__ == '__main__':
    unittest.main()



if __name__ == '__main__':
    unittest.main()