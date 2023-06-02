"""

File:    ImportTests.py

     Some test cases ..

"""
import unittest

from wwpdb.utils.emdb.EmdbSchema import EmdbSchema  # noqa: F401 pylint: disable=unused-import
from wwpdb.utils.emdb.cif_emdb_translator.cif_emdb_translator import CifEMDBTranslator  # noqa: E402
from wwpdb.utils.emdb.atomcheck.atomcheck import read_args  # noqa: F401 pylint: disable=unused-import


class ImportTests(unittest.TestCase):
    def setUp(self):
        pass

    def testInstantiate(self):
        _ce = CifEMDBTranslator()  # noqa: F841
        _es = EmdbSchema()  # noqa: F841


if __name__ == "__main__":
    unittest.main()
