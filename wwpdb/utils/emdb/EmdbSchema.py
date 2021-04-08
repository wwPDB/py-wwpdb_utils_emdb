##
# File:    EmdbSchema.py
# Date:    29-March-2021
"""
Wrapper for access to EMDB schema.

The schema is contained in this package
"""

import os


class EmdbSchema(object):
    def getSchemaPath(self):
        spath = os.path.join(os.path.dirname(__file__), "data", "emdb-v3.xsd")
        return os.path.abspath(spath)
