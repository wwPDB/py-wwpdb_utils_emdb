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

import platform
import os
import unittest

#####################  setup DepUi test environment ############
HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
TESTOUTPUT = os.path.join(HERE, 'test-output', platform.python_version())
if not os.path.exists(TESTOUTPUT):
    os.makedirs(TESTOUTPUT)
mockTopPath = os.path.join(TOPDIR, 'wwpdb', 'mock-data')
rwMockTopPath = os.path.join(TESTOUTPUT)

# Must create config file before importing ConfigInfo
from wwpdb.utils.testing.SiteConfigSetup  import SiteConfigSetup
from wwpdb.utils.testing.CreateRWTree import CreateRWTree
# Copy site-config and selected items
crw = CreateRWTree(mockTopPath, TESTOUTPUT)
crw.createtree(['site-config', 'depuiresources'])
# Use populate r/w site-config using top mock site-config
SiteConfigSetup().setupEnvironment(rwMockTopPath, rwMockTopPath)

# Setup DepUI specific directories
from wwpdb.utils.config.ConfigInfo import ConfigInfo
import os
import os.path
cI = ConfigInfo()
FILE_UPLOAD_TEMP_DIR = os.path.join(
    cI.get("SITE_DEPOSIT_STORAGE_PATH"),
    "deposit",
    "temp_files")
if not os.path.exists(FILE_UPLOAD_TEMP_DIR):
    os.makedirs(FILE_UPLOAD_TEMP_DIR)

# Django envivonment setup
os.environ['DJANGO_SETTINGS_MODULE'] = "wwpdb.apps.deposit.settings"
os.environ['IN_ANNOTATION'] = "no"
##################################################

from wwpdb.utils.emdb.cif_emdb_translator.cif_emdb_translator import CifEMDBTranslator


class ImportTests(unittest.TestCase):
    def setUp(self):
        pass

    def testInstantiate(self):
        cT = CifEMDBTranslator()
        pass

    
