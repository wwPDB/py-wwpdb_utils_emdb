##
# File: ImportTests.py
# Date:  06-Oct-2018  E. Peisach
#
# Updates:
##
"""Test cases for emdb - simply import everything to ensure imports work"""

__docformat__ = "restructuredtext en"
__author__ = "Ezra Peisach"
__email__ = "peisach@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import unittest

from wwpdb.utils.emdb.cif_emdb_translator.cif_emdb_translator import CifEMDBTranslator

class ImportTests(unittest.TestCase):
    def setUp(self):
        pass

    def testInstantiate(self):
        cT = CifEMDBTranslator()
        pass

    
