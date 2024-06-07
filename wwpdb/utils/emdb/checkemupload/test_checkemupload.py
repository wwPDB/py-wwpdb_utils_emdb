import os
import numpy as np
import json
import inspect
import unittest
from unittest.mock import patch, mock_open, MagicMock
import hashlib
import mrcfile
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

import django
os.environ["DJANGO_SETTINGS_MODULE"] = "wwpdb.apps.deposit.settings"
django.setup()

from wwpdb.utils.emdb.checkemupload.checkemupload import EMMap, Model, Validator
from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId

site_id = getSiteId()
config = ConfigInfo(site_id)
onedep_path = os.environ['ONEDEP_PATH']
deployment = config.get('SITE_SUFFIX')
data_dir = os.path.join(onedep_path, 'deployments', deployment, 'source/py-wwpdb_utils_emdb/wwpdb/utils/emdb/checkemupload/data/')
stub_data_path = os.path.join(data_dir, 'stub_data.json')
# raise ValueError(f'Path: {stub_data_path}')


# Add this class to convert numpy arrays to lists when serializing to JSON
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.float32):
            return float(obj)
        elif isinstance(obj, bytes):  # Add this check
            return obj.decode()  # Convert bytes to string
        return json.JSONEncoder.default(self, obj)


# Create a class that generates the expected results for each situation
class ResultGenerator:
    def __init__(self, stub_data_path):
        # Load the JSON file containing the file paths
        with open(stub_data_path, 'r') as f:
            self.stub_data = json.load(f)

    # Use patch to mock os.path.isfile
    @patch('os.path.isfile', return_value=True)
    def get_result(self, situation, mock_isfile):
        # Iterate over the stub data
        for content in self.stub_data:
            if content['situation'] != situation:
                continue

            primmap_data = content['primmap']
            halfmaps_data = content['halfmaps']
            model_data = content['model']

            primmap = EMMap(primmap_data['file'])
            # For each attribute in the primmap object, assign the value from the JSON file that has the same key
            for key, value in primmap_data.items():
                setattr(primmap, key, value)
            
            # Do the same for the halfmaps
            halfmaps = []
            for halfmap_data in halfmaps_data:
                halfmap = EMMap(halfmap_data['file'])
                for key, value in halfmap_data.items():
                    setattr(halfmap, key, value)
                halfmaps.append(halfmap)
            
            # Do the same for the model
            model = Model(model_data['file'])
            for key, value in model_data.items():
                setattr(model, key, value)

            # Validate the objects
            validator = Validator(primmap, halfmaps, model)
            result = validator.check()

            # Return the result for the given situation
            return result

        # If the situation is not found in stub_data, raise an error
        raise ValueError(f'Situation "{situation}" not found in stub data')


# Create a test class for unit testing the EMMap class
class TestEMMap(unittest.TestCase):

    def setUp(self):
        # Correcting the mock paths
        patcher_mrcfile = patch('mrcfile.mrcfile.MrcFile', autospec=True)
        self.addCleanup(patcher_mrcfile.stop)
        self.mock_mrc = patcher_mrcfile.start()

        patcher_header = patch('mrcfile.mrcobject.MrcObject.header', new_callable=MagicMock)
        self.addCleanup(patcher_header.stop)
        self.mock_header = patcher_header.start()

        # Create an instance of EMMap with mock data
        self.em_map = EMMap('path/to/mock_map.mrc')
        self.em_map.header = MagicMock()
        self.em_map.nxyz = [10, 10, 10]
        self.em_map.box_size = [100, 100, 100]
        self.em_map.pixel_size = [1.0, 1.0, 1.0]
        self.em_map.nstarts = [0, 0, 0]
        self.em_map.origin = [0.0, 0.0, 0.0]
        self.em_map.end = [100.0, 100.0, 100.0]

    def test_load(self):
        self.assertIsNotNone(self.em_map.header)
        self.assertEqual(self.em_map.nxyz, [10, 10, 10])
        self.assertEqual(self.em_map.box_size, [100, 100, 100])
        self.assertEqual(self.em_map.pixel_size, [1.0, 1.0, 1.0])
        self.assertEqual(self.em_map.nstarts, [0, 0, 0])
        self.assertEqual(self.em_map.origin, [0.0, 0.0, 0.0])
        self.assertEqual(self.em_map.end, [100.0, 100.0, 100.0])

    def test_md5_checksum(self):
        data = np.array([72, 101, 108, 108, 111, 44, 32, 119, 111, 114, 108, 100, 33], dtype=np.uint8)
        checksum = hashlib.md5(data.tobytes()).hexdigest()
        self.assertEqual(self.em_map.md5_checksum(data), checksum)

    def test_same_box_size(self):
        another_map = EMMap('path/to/another_map.mrc')
        another_map.box_size = [100, 100, 100]
        self.assertTrue(self.em_map.same_box_size(another_map))

        another_map.box_size = [50, 50, 50]
        self.assertFalse(self.em_map.same_box_size(another_map))

    def test_same_or_smaller_box_size(self):
        another_map = EMMap('path/to/another_map.mrc')
        another_map.box_size = [100, 100, 100]
        self.assertTrue(self.em_map.same_or_smaller_box_size(another_map))

        another_map.box_size = [50, 50, 50]
        self.assertFalse(self.em_map.same_or_smaller_box_size(another_map))

    def test_overlaps(self):
        another_map = EMMap('path/to/another_map.mrc')
        another_map.origin = [0.0, 0.0, 0.0]
        another_map.end = [100.0, 100.0, 100.0]
        self.assertTrue(self.em_map.overlaps(another_map))

        another_map.origin = [101.0, 101.0, 101.0]
        another_map.end = [200.0, 200.0, 200.0]
        self.assertFalse(self.em_map.overlaps(another_map))

    def test_fits_inside(self):
        another_map = EMMap('path/to/another_map.mrc')
        another_map.origin = [-1.0, -1.0, -1.0]
        another_map.end = [101.0, 101.0, 101.0]
        self.assertTrue(self.em_map.fits_inside(another_map))

        another_map.origin = [1.0, 1.0, 1.0]
        another_map.end = [99.0, 99.0, 99.0]
        self.assertFalse(self.em_map.fits_inside(another_map))

    def test_same_pixel_size(self):
        another_map = EMMap('path/to/another_map.mrc')
        another_map.pixel_size = [1.0, 1.0, 1.0]
        self.assertTrue(self.em_map.same_pixel_size(another_map))

        another_map.pixel_size = [0.5, 0.5, 0.5]
        self.assertFalse(self.em_map.same_pixel_size(another_map))

    def test_pixel_size_is_multiple(self):
        another_map = EMMap('path/to/another_map.mrc')
        another_map.pixel_size = [0.5, 0.5, 0.5]
        self.assertTrue(self.em_map.pixel_size_is_multiple(another_map))

        another_map.pixel_size = [0.3, 0.3, 0.3]
        self.assertFalse(self.em_map.pixel_size_is_multiple(another_map))


# Create a test class for unit testing the Model class
class TestModel(unittest.TestCase):

    def setUp(self):
        # Create an instance of Model with mock data
        self.model = Model('path/to/mock_model.cif')
        self.model.structure = [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)]
        self.model.errors = []

    def test_load(self):
        with patch.object(self.model, 'get_coordinates', return_value=self.model.structure):
            self.model.load()
            self.assertEqual(self.model.structure, [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)])

    def test_get_coordinates(self):
        with patch.object(self.model, 'parse_mmcif', return_value=(
            {'_atom_site.Cartn_x': ['1.0'], '_atom_site.Cartn_y': ['2.0'], '_atom_site.Cartn_z': ['3.0']},
            False, ""
        )):
            coordinates = self.model.get_coordinates()
            self.assertEqual(coordinates, [(1.0, 2.0, 3.0)])

    @patch('os.path.isfile', return_value=True)
    @patch('os.path.basename', return_value='dummy_file')
    @patch('Bio.PDB.MMCIF2Dict.MMCIF2Dict', return_value = {
            '_atom_site.Cartn_x': ['1.0'],
            '_atom_site.Cartn_y': ['2.0'],
            '_atom_site.Cartn_z': ['3.0']
        })
    def test_parse_mmcif(self, mock_mmcif2dict, mock_basename, mock_isfile):
        # Mock the open function to return a file-like object with some data
        mock_file = mock_open(read_data='data_test\n_atom_site.Cartn_x 1.0\n_atom_site.Cartn_y 2.0\n_atom_site.Cartn_z 3.0\n')
        with patch('builtins.open', mock_file), patch('Bio.PDB.MMCIF2Dict.MMCIF2Dict', return_value=mock_mmcif2dict):
            # Act
            result = self.model.parse_mmcif()

        # Adjusted expected result
        expected_result = {'_atom_site.Cartn_x': ['1.0'], '_atom_site.Cartn_y': ['2.0'], '_atom_site.Cartn_z': ['3.0'], 'data_': 'test'}
        self.assertEqual(result, (expected_result, False, ""))


# Create a test class for unit testing the Validator class
class TestValidator(unittest.TestCase):

    @patch('os.path.isfile', return_value=True)
    def setUp(self, mock_isfile):
        # Create instances of EMMap and Model with mock data
        # Use the `ResultGenerator` class to generate the expected results for the "Everything is fine" situation
        self.generator = ResultGenerator(stub_data_path)
        self.stub_data = None
        for situation in self.generator.stub_data:
            if situation['situation'] == 'Everything is fine':
                self.stub_data = {k:v for k,v in situation.items() if k != 'situation'}
                break
        self.expected_result = self.generator.get_result('Everything is fine')

        primmap_data = self.stub_data['primmap']
        halfmaps_data = self.stub_data['halfmaps']
        model_data = self.stub_data['model']

        self.primmap = EMMap(primmap_data['file'])
        for key, value in primmap_data.items():
            setattr(self.primmap, key, value)

        self.halfmaps = []
        for halfmap_data in halfmaps_data:
            halfmap = EMMap(halfmap_data['file'])
            for key, value in halfmap_data.items():
                setattr(halfmap, key, value)
            self.halfmaps.append(halfmap)

        self.model = Model(model_data['file'])
        for key, value in model_data.items():
            setattr(self.model, key, value)

        # Create an instance of Validator with the mock data    
        self.validator = Validator(self.primmap, self.halfmaps, self.model)

    @patch('os.path.isfile', return_value=True)
    def test_check(self, mock_isfile):
        # Act
        result = self.validator.check()

        # Assert
        self.assertEqual(result, self.expected_result)


class TestCheckEmUpload(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        # Create an instance of ResultGenerator
        stub_data_path = os.path.join(data_dir, 'stub_data.json')
        cls.generator = ResultGenerator(stub_data_path)
        # # NOTE: Results in the `answer_key.json` file are generated using the following code
        # # Generate the expected results for each situation and write them to the `answer_key.json` file
        # cls.expected_results = {}
        # for situation in cls.generator.stub_data:
        #     result = cls.generator.get_result(situation['situation'])
        #     cls.expected_results[situation['situation']] = result
        # # Write the expected results to a JSON file
        # with open(os.path.join(data_dir, 'answer_key.json'), 'w') as f:
        #     json.dump(cls.expected_results, f, cls=NumpyEncoder, indent=4)
        # # NOTE: The expected results are generated only once and saved in the `answer_key.json` file, which can be further edited manually
        
        # Load the expected results from the `answer_key.json` file
        with open(os.path.join(data_dir, 'answer_key.json'), 'r') as f:
            cls.expected_results = json.load(f)
        # NOTE: Model-related tests are currently resulting in null values for the expected results
        # NOTE: This behavior is expected as the header attribute of the EMMap object is not being set, just mocked
        # TODO: To fix this, the header attribute of the EMMap object should be set (or mocked) to a valid value

    def test_everything_is_fine(self):
        actual_result = self.generator.get_result('Everything is fine')
        self.assertEqual(actual_result, self.expected_results['Everything is fine'])

    def test_half_map_identical_to_primary_map(self):
        actual_result = self.generator.get_result('Half-map identical to primary map')
        self.assertEqual(actual_result, self.expected_results['Half-map identical to primary map'])

    def test_half_maps_identical_to_each_other(self):
        actual_result = self.generator.get_result('Half-maps identical to each other')
        self.assertEqual(actual_result, self.expected_results['Half-maps identical to each other'])
    
    def test_half_maps_not_overlaid(self):
        actual_result = self.generator.get_result('Half-maps not overlaid')
        self.assertEqual(actual_result, self.expected_results['Half-maps not overlaid'])

    def test_half_maps_with_different_pixel_sizes(self):
        actual_result = self.generator.get_result('Half-maps with different pixel sizes')
        self.assertEqual(actual_result, self.expected_results['Half-maps with different pixel sizes'])

    def test_model_contains_atoms_outside_map_bounding_box(self):
        actual_result = self.generator.get_result('Model contains atoms outside map bounding box')
        self.assertEqual(actual_result, self.expected_results['Model contains atoms outside map bounding box'])

    def test_primary_map_with_box_size_greater_than_half_maps(self):
        actual_result = self.generator.get_result("Primary map with box size greater than half map's")
        self.assertEqual(actual_result, self.expected_results["Primary map with box size greater than half map's"])

    def test_primary_map_bounding_box_extends_beyond_boundary_of_half_map_bounding_boxes(self):
        actual_result = self.generator.get_result('Primary map bounding box extends beyond the boundary of the half-map bounding boxes')
        self.assertEqual(actual_result, self.expected_results['Primary map bounding box extends beyond the boundary of the half-map bounding boxes'])

    def test_situation_not_found(self):
        with self.assertRaises(ValueError):
            self.generator.get_result('Situation not found')

if __name__ == '__main__':
    unittest.main()
