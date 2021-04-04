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

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
TESTOUTPUT = os.path.join(HERE, 'test-output', platform.python_version())
if not os.path.exists(TESTOUTPUT):
    os.makedirs(TESTOUTPUT)
mockTopPath = os.path.join(TOPDIR, 'wwpdb', 'mock-data')

# Must create config file before importing ConfigInfo
from wwpdb.utils.testing.SiteConfigSetup  import SiteConfigSetup
SiteConfigSetup().setupEnvironment(TESTOUTPUT, mockTopPath)

#
from wwpdb.utils.emdb.cif_emdb_translator.cif_emdb_translator import CifEMDBTranslator


class ImportTests(unittest.TestCase):
    def setUp(self):
        self.__inpfile = os.path.join(mockTopPath, 'EMD', 'emd-0000.cif')
        self.__outfile = os.path.join(TESTOUTPUT, 'emd-0000.xml')
        self.__logfile = os.path.join(TESTOUTPUT, 'emd-0000.log')
        pass

    @staticmethod
    def testInstantiate():
        """Tests simple instantiation"""
        CifEMDBTranslator()

    def testTranslateSuppressed(self):
        """Tests translation of suppressed input"""

        # Changed to not specify schema - let system determine
        translator = CifEMDBTranslator()
        translator.set_logger_logging(log_error=True, error_log_file_name=self.__logfile)
        translator.read_emd_map_v2_cif_file()
        translator.translate_and_validate(in_cif=self.__inpfile, out_xml=self.__outfile)
        # This will close the output file
        translator.write_logger_logs(write_error_log=True)

        self.assertTrue(translator.is_translation_log_empty, 'Translator failed')
        self.assertTrue(os.path.exists(self.__outfile), "No output file")

if __name__ == "__main__":
    unittest.main()
