#!/usr/bin/env python
"""
unittest.py

Test that the translator is working
"""


__author__ = "Ardan Patwardhan, Sanja Abbott"
__email__ = "ardan@ebi.ac.uk, sanja@ebi.ac.uk"
__date__ = "2018-06-28"

import os.path
import unittest
import glob
# from wwpdb.utils.emdb.cif_emdb_translator.cif_emdb_translator import CifEMDBTranslator
from wwpdb.utils.emdb.cif_emdb_translator.cif_emdb_legacy_translator import CifEMDBTranslator


class TestTranslator(unittest.TestCase):
    """Run the translator on a bunch of example files"""

    testOutputXmlDir = "data/test/xml_v3_out"
    testInputCifDir = "data/cif"
    schema = "emdb.xsd"

    def test_cif2xml(self):

        # create the translator object
        translator = CifEMDBTranslator()
        translator.set_logger_logging(True, True, True, False)

        test_cifs = glob.glob(self.testInputCifDir + "/*.cif")
        for test_cif in test_cifs:
            print("In file: " + test_cif)
            out_filename = os.path.splitext(os.path.basename(test_cif))[0]
            xml_out = os.path.join(self.testOutputXmlDir + "/" + out_filename + ".xml")
            translator.translate_and_validate(test_cif, xml_out, self.schema)
            translator.write_logger_logs(True, True, True)


if __name__ == "__main__":
    unittest.main()
