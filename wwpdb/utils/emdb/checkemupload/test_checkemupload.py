import json
from wwpdb.utils.emdb.checkemupload.checkemupload import EMMap, Model
import inspect

path = '/nfs/production/gerard/pdbe/onedep/deployments/emdb_dev_1/source/py-wwpdb_utils_emdb/wwpdb/utils/emdb/checkemupload/data/stub_objects.json'

# Load the JSON file containing the file paths
with open(path, 'r') as f:
    data = json.load(f)

# Iterate over the data
for content in data:
    situation = content['situation']
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

    # Print the values of the attributes of the objects
    print(f'Situation: {situation}')
    print(f'Primmap: {json.dumps({k: v for k, v in primmap.__dict__.items() if not inspect.isfunction(v)}, indent=4)}')
    print(f'Halfmaps: {json.dumps([{k: v for k, v in halfmap.__dict__.items() if not inspect.isfunction(v)} for halfmap in halfmaps], indent=4)}')
    print(f'Model: {json.dumps({k: v for k, v in model.__dict__.items() if not inspect.isfunction(v)}, indent=4)}')
    print('\n\n')