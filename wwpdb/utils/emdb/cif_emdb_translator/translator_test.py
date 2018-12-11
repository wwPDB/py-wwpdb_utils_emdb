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
__date__ = '2017-08-24'

import os.path
import logging
import unittest

from lxml import etree
import csv
import sys
import glob

from getLatest import PDBprocessedWhere
from em_emd_conversion import convert
from cif_emdb_translator import CifEMDBTranslator

class TestTranslator(unittest.TestCase):
    """Run the translator on a bunch of example files"""
    pdbeInputCifDir = 'data/cif'
    rcsbInputCifDir = 'data/rcsb_cifs'
    outputXmlDir = 'data/test/xml_v3_out'
    testInputCifDir = 'data/test/cif'

    schema = "emdb.xsd"

    process_rcsb = False
    process_pdbe = True
    process_test = False

    def test_cif2xml(self):

        # create the translator object
        translator = CifEMDBTranslator()
        # translator.set_logger_logging(True, True, True, False)
        translator.set_logger_logging(log_error=True, error_log_file_name='/nfs/msd/work2/sanja/cif_emdb_translator/ERROR.log')
        # translator.set_show_log_id(True)
        # Reads mmcif_pdbx_v5_next.dic that contains information
        # about how the em categories map to the emd categories
        #translator.read_emd_map_v2_cif_file()
        i = 1

        if self.process_test:
            try:
                test_cifs = glob.glob(self.testInputCifDir+'/*.cif')
                for test_cif in test_cifs:
                    if test_cif.find('D_1000232975_emdb.cif') != -1:
                        xml_out = os.path.join(self.testInputCifDir+'/test.xml')
                        translator.translate_and_validate(test_cif, xml_out, self.schema)
                        translator.write_logger_logs(True, True, True)
                        print
                        print 'TEST CLASS LOGS'
                        print
                        if translator.is_translation_log_empty:
                            print 'NO ERRORS FOUND'
                        else:
                            for entry_log in translator.translation_log.logs:
                                print entry_log.id
                                if entry_log.is_error_log_empty:
                                    print
                                    print "NO ERRORS FOUND FOR THIS ENTRY"
                                    print
                                else:
                                    for err in entry_log.errors:
                                        print
                                        print err.log_text
                                        print
            except Exception as ex:
                print ex

                #     if self.process_rcsb:
                #         try:
                #             # lists of EMD ids and their corresponding D_ ids - taken from the database dump
                #             #listOfEMIDs,  listOfDepIDs = self.getAllIDs()
                #             # list of actual files
                #             rcsb_cifs = glob.glob(self.rcsbInputCifDir+'/*.cif')
                #             for rcsb_cif in rcsb_cifs:
                #                 print i
                #                 i = i + 1
                #                 if rcsb_cif.find('D_') != -1:
                #                     # find EMD id for this file
                #                     d_index = rcsb_cif.index('D_')
                #                     d_id = rcsb_cif[d_index : d_index+12]
                #                     xml_out = os.path.join(self.outputXmlDir, d_id + '.xml')
                #                     translator.translate_and_validate(rcsb_cif, xml_out, self.schema)
                #         except Exception as ex:
                #             print ex
                #
        if self.process_pdbe:
            try:
                listOfEMIDs,  listOfDepIDs = self.getAllIDs()
                if listOfEMIDs:
                    print len(listOfEMIDs)
                    j = 0
                    for id in listOfEMIDs:
                        depID = listOfDepIDs[j]
                        j = j + 1
                        #print depID
                        #if depID == 'D_1200000799' or depID == 'D_1200005141':
                        if id == "EMD-8142":# or id == "EMD-8057":
                            print id
                            try:
                                fileType, f, copyStatic = PDBprocessedWhere(id).Extension()
                                if f:
                                    print i
                                    i = i + 1
                                    #Check if the cif file is in the _emd space
                                    if not self.isCifInEMDSpace(f):
                                        conv_f = os.path.join(self.pdbeInputCifDir, id + '.cif')
                                        print conv_f
                                        # Check if the converted file exits
                                        if not self.convertedCifToEMDSpaceExits(conv_f):
                                            # There is no converted file, convert it now
                                            convert(f, conv_f).em2emd()
                                    if os.path.exists(conv_f) and os.path.isfile(conv_f):
                                        of = os.path.join(self.outputXmlDir, id + '.xml')
                                        print of
                                        print "conv_f %s" % conv_f
                                        translator.translate_and_validate(conv_f, of, self.schema)
                                        a_log = translator.current_entry_log
                                        log_id = a_log.id
                                    else:
                                        print 'The file ' + conv_f + 'cannot be converted into the _emd space and therefore, cannot be translated'
                            except IOError as exp:
                                print exp
                            print
                            print 'PDBE LOGGER LOGS'
                            print
                            translator.write_logger_logs(write_error_log=True)
                print
                print 'PDBE CLASS LOGS'
                print
                if translator.is_translation_log_empty:
                    print 'NO ERRORS FOUND'
                else:
                    for entry_log in translator.translation_log.logs:
                        print entry_log.id
                        if entry_log.is_error_log_empty:
                            print
                            print "NO ERRORS FOUND FOR THIS ENTRY"
                            print
                        else:
                            for err in entry_log.errors:
                                print
                                print err.log_text
                                print
            except IOError as exp:
                print exp

    def isCifInEMDSpace(self, f):
        if '_emd' in open(f).read():
            return True
        else:
            return False

    def convertedCifToEMDSpaceExits(self, f):
        if os.path.exists(f) and os.path.isfile(f):
            return True
        else:
            return False

    def getAllIDs(self):
        ids = []
        d_ids = []
        latestEMDBDump = '/nfs/pdbe_da/production/data/for_release/emd/em_db_status.csv'
        with open(latestEMDBDump) as csvfile:
            reader = csv.DictReader(csvfile)
            for col in reader:
                ids.append(col['emdb_id'])
                d_ids.append(col['dep_id'])
        return ids, d_ids

    def validateWithSchemaV20(self,headerToValidate):

        try:
            f = open(self.schema, 'r')
        except:
            return False
        else:
            schema_root = etree.XML(f.read())
            theSchema = etree.XMLSchema(schema_root)
            xml_parser = etree.XMLParser(schema=theSchema)
            validate = self.validateFile(xml_parser, headerToValidate)

            if validate:
                #print "%s validates" % headerToValidate
                return True
            else:
                #print "%s FAILED validation" % headerToValidate
                return False
        finally:
            f.close()

    def validateFile(self, theParser, xmlfilename):
        #this method can be used to validate any schema against any file
        try:
            f = open(xmlfilename, 'r')
            try:
                etree.fromstring(f.read(), theParser)
            except etree.XMLSyntaxError as err:
                #print "Can't parse %s file. Error is %s" % (xmlfilename, err)
                return False
            except etree.XMLSchemaError as err:
                #print "XMLSchemaError: %s" % err
                return False
            return True
        except:
            #print "Can't open: %s file" % xmlfilename
            return False
        finally:
            f.close()

if __name__ == '__main__':
    unittest.main()
