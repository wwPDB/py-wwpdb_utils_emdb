import os
import django
os.environ["DJANGO_SETTINGS_MODULE"] = "wwpdb.apps.deposit.settings"
django.setup()

import simplejson

# from wwpdb.apps.deposit.depui.FileConversion import FileConversion
from wwpdb.apps.deposit.depui.common_functions import getDir, getFile
# from wwpdb.utils.wf.dbapi.WfDbApi import WfDbApi

# username = "D_800299"
# output_dir = getDir(username, "tempdep")
# print(f"\noutput_dir: {output_dir}")

# wf_op = "emmapcheck"
# fc = FileConversion()
# wfApi = WfDbApi(verbose=True)
# status = fc.wfStart(username, wfApi, wf_op, "wf_op_emmapcheck_fs_deposit.xml")
# wf_data = fc.wfStatus(username, wfApi)
# finished = fc.wfFinished(username, wfApi, wf_data, wf_op)

# if status == 0 or not finished:
#     print("Unable to perform checks on uploaded files!")

# path2output = getFile(username, "tempdep", "em-map-report", "json", versionId='latest', mileStone=None, wfInstanceId=None, partNumber=None)

# try:
#     with open(path2output, 'r', encoding='utf-8') as f:
#         output = simplejson.load(f)
# except Exception as e:
#     print(f"\nException: {e}")

# # Iterate over the files in the output directory sorted by descending modification time
# for root, dirs, files in sorted(os.walk(output_dir), key=lambda x: os.stat(os.path.join(x[0], x[2][0])).st_mtime, reverse=True):
#     level = root.replace(output_dir, '').count(os.sep)
#     indent = ' ' * 4 * (level)
#     print(f"{indent}{os.path.basename(root)}/")
#     subindent = ' ' * 4 * (level + 1)
#     for f in files:
#         print(f"{subindent}{f}")


import unittest
from unittest.mock import patch, MagicMock

class TestFileConversion(unittest.TestCase):
    @patch('os.stat')
    @patch('os.walk')
    @patch('wwpdb.apps.deposit.depui.FileConversion.FileConversion')
    @patch('wwpdb.utils.wf.dbapi.WfDbApi.WfDbApi')
    def test_file_conversion(self, mock_WfDbApi, mock_FileConversion, mock_walk, mock_stat):
        os.environ["DJANGO_SETTINGS_MODULE"] = "wwpdb.apps.deposit.settings"
        django.setup()

        username = "D_800299"
        output_dir = getDir(username, "tempdep")

        wf_op = "emmapcheck"
        fc = mock_FileConversion.return_value
        wfApi = mock_WfDbApi.return_value
        status = fc.wfStart(username, wfApi, wf_op, "wf_op_emmapcheck_fs_deposit.xml")
        wf_data = fc.wfStatus(username, wfApi)
        finished = fc.wfFinished(username, wfApi, wf_data, wf_op)

        self.assertNotEqual(status, 0)
        self.assertTrue(finished)

        path2output = getFile(username, "tempdep", "em-map-report", "json", versionId='latest', mileStone=None, wfInstanceId=None, partNumber=None)

        with patch('builtins.open', new=MagicMock()) as mock_open:
            mock_open.return_value.__enter__.return_value = MagicMock()
            mock_open.return_value.__enter__.return_value.read.return_value = '{"key": "value"}'
            with open(path2output, 'r', encoding='utf-8') as f:
                output = simplejson.load(f)
            self.assertEqual(output, {"key": "value"})

if __name__ == '__main__':
    unittest.main()