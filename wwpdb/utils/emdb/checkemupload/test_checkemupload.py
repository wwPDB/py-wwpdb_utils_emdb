import os
import numpy as np
import json
import inspect
import unittest
from unittest.mock import patch

import django
os.environ["DJANGO_SETTINGS_MODULE"] = "wwpdb.apps.deposit.settings"
django.setup()

from wwpdb.utils.emdb.checkemupload.checkemupload import EMMap, Model, Validator
from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId

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

site_id = getSiteId()
config = ConfigInfo(site_id)
onedep_path = os.environ['ONEDEP_PATH']
deployment = config.get('SITE_SUFFIX')
data_dir = os.path.join(onedep_path, 'deployments', deployment, 'source/py-wwpdb_utils_emdb/wwpdb/utils/emdb/checkemupload/data/')
stub_data_path = os.path.join(data_dir, 'stub_data.json')
# raise ValueError(f'Path: {stub_data_path}')

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


class TestCheckEmUpload(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Load the expected results from answer_key.json
        results_path = os.path.join(data_dir, 'answer_key.json')
        with open(results_path, 'r') as f:
            cls.expected_results = json.load(f)

        # Create an instance of ResultGenerator
        stub_data_path = os.path.join(data_dir, 'stub_data.json')
        cls.generator = ResultGenerator(stub_data_path)

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
