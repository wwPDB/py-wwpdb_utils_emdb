#!/usr/bin/env python
"""
unittest.py

Test that the translator is working

When the input cifs are loaded from /wwpdb_da/da_top/data/archive the processing must be done on triton-2!!!!
The cifs from the archive are in _em space so they need to be converted into the _emd space.

TODO:
1) Should discover files in data/cif


Version history:
0.1, 2015-07-26, Ardan Patwardhan: Test on files in data/cif


Copyright [2015] EMBL - European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the
"License"); you may not use this file except in
compliance with the License. You may obtain a copy of
the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on
an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
"""


__author__ = 'Ardan Patwardhan, Sanja Abbott'
__email__ = 'ardan@ebi.ac.uk, sanja@ebi.ac.uk'
__date__ = '2018-06-28'

import os.path
import unittest
import glob
from cif_emdb_translator import CifEMDBTranslator


class TestTranslator(unittest.TestCase):
    """Run the translator on a bunch of example files"""
    testOutputXmlDir = 'data/test/xml_out'
    testInputCifDir = 'data/test/cif'
    schema = "emdb.xsd"

    def test_cif2xml(self):

        # create the translator object
        translator = CifEMDBTranslator()
        translator.set_logger_logging(True, True, True, False)
        i = 1

        test_cifs = glob.glob(self.testInputCifDir+'/*.cif')
        for test_cif in test_cifs:
            out_filename = os.path.splitext(os.path.basename(test_cif))[0]
            xml_out = os.path.join(self.testOutputXmlDir+'/'+out_filename+'.xml')
            translator.translate_and_validate(test_cif, xml_out, self.schema)
            translator.write_logger_logs(True, True, True)


if __name__ == '__main__':
    unittest.main()
