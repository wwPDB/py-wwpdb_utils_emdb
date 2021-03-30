##
# File: SchemaPathTests.py
# Date:  06-Oct-2018  E. Peisach
#
# Updates:
##
"""Test cases for emdb schema path"""

import os
import unittest

from wwpdb.utils.emdb.EmdbSchema import EmdbSchema


class SchamePathTests(unittest.TestCase):
    def testSchemaPath(self):
        es = EmdbSchema()
        spath = es.getSchemaPath()

        ok = os.path.exists(spath)
        self.assertTrue(ok, "EMDB schema missing: %s" % spath)
        
