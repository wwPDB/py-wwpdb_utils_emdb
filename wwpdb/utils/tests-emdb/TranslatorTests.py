##
# File: TranslatorTests.py
# Date:  06-Oct-2018  E. Peisach
#
# Updates:
##
"""Test cases for emdb - simply import everything to ensure imports work"""

__docformat__ = "restructuredtext en"
__author__ = "Ezra Peisach"
__email__ = "peisach@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"

import platform
import os
import unittest

try:
    from unittest.mock import MagicMock, patch
except ImportError:
    from mock import MagicMock, patch

from wwpdb.utils.emdb.cif_emdb_translator.cif_emdb_translator import CifEMDBTranslator  # noqa: E402


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
TESTOUTPUT = os.path.join(HERE, "test-output", platform.python_version())
if not os.path.exists(TESTOUTPUT):
    os.makedirs(TESTOUTPUT)


class ImportTests(unittest.TestCase):
    def setUp(self):
        self.__inpfile = os.path.join(HERE, "data", "emd-0000.cif")
        self.__outfile = os.path.join(TESTOUTPUT, "emd-0000.xml")
        self.__logfile = os.path.join(TESTOUTPUT, "emd-0000.log")

    @staticmethod
    def testInstantiate():
        """Tests simple instantiation"""
        CifEMDBTranslator()

    @patch("wwpdb.utils.emdb.cif_emdb_translator.cif_emdb_translator.ConfigInfoAppEm")
    def testTranslateSuppressed(self, ciamock):
        """Tests translation of suppressed input"""
        instance = MagicMock()
        instance.get_emd_mapping_file_path.return_value = os.path.join(HERE, "data", "emd", "emd_map_v2.cif")
        ciamock.return_value = instance

        # Changed to not specify schema - let system determine
        translator = CifEMDBTranslator()
        translator.set_logger_logging(log_error=True, error_log_file_name=self.__logfile)
        translator.read_emd_map_v2_cif_file()
        translator.translate_and_validate(in_cif=self.__inpfile, out_xml=self.__outfile)
        # This will close the output file
        translator.write_logger_logs(write_error_log=True)

        self.assertTrue(translator.is_translation_log_empty, "Translator failed")
        self.assertTrue(os.path.exists(self.__outfile), "No output file")


if __name__ == "__main__":
    unittest.main()
