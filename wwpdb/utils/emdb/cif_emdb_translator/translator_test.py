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
__date__ = '2019-05-17'

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
    testInputCifDir = 'data/cif'
    db_dump = '/nfs/msd/em/dep_id_to_emdb_id.csv'

    schema = "emdb.xsd"

    process_rcsb = False
    process_pdbe = True
    process_test = False

    def test_cif2xml(self):

        # create the translator object
        translator = CifEMDBTranslator()
        # get the script's location
        script_loc = os.path.dirname(os.path.realpath(__file__))
        # translator.set_logger_logging(True, True, True, False)
        error_log_loc = os.path.join(script_loc, 'ERROR.log')
        translator.set_logger_logging(log_error=True, error_log_file_name=error_log_loc)
        # translator.set_show_log_id(True)
        # Reads mmcif_pdbx_v5_next.dic that contains information
        # about how the em categories map to the emd categories
        #translator.read_emd_map_v2_cif_file()
        i = 1

        if self.process_test:
            try:
                test_dir = os.path.join(script_loc, self.testInputCifDir)
                test_cifs = glob.glob(test_dir + '/*.cif')
                for test_cif in test_cifs:
                    emd_id = os.path.splitext(os.path.basename(test_cif))[0]
                    out_name = emd_id + '.xml'
                    if True: # emd_id == 'EMD-10119':
                        print 'translating %s' % test_cif
                        xml_out = os.path.join(self.outputXmlDir, out_name)
                        translator.translate_and_validate(test_cif, xml_out, self.schema)
                        translator.write_logger_logs(True, True, True)
                        #print '\nTEST CLASS LOGS\n'
                        if translator.is_translation_log_empty:
                            print '\nNO ERRORS FOUND\n'
                        else:
                            for entry_log in translator.translation_log.logs:
                                print 'entry_log.id %s' % entry_log.id
                                if entry_log.is_error_log_empty:
                                    print "\nNO ERRORS FOUND FOR THIS ENTRY\n"
                                else:
                                    for err in entry_log.errors:
                                        print '\nERRORS: %s' % err.log_text
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
                    print "# of EM IDs - %s\n" % len(listOfEMIDs)
                    j = 0
                    for id in listOfEMIDs:
                        depID = listOfDepIDs[j]
                        j = j + 1
                        #print "dep_id: %s" % depID
                        #if depID == 'D_1200000799' or depID == 'D_1200005141':
                        if True: #id == "EMD-10119":# or id == "EMD-8057":
                            # print "EMD ID: %s" % id
                            try:
                                fileType, f, copyStatic = PDBprocessedWhere(id).Extension()
                                if f:
                                    print '%s: %s' % (i, id)
                                    i = i + 1
                                    #Check if the cif file is in the _emd space
                                    if not self.isCifInEMDSpace(f):
                                        conv_f = os.path.join(self.pdbeInputCifDir, id + '.cif')
                                        #print conv_f
                                        # Check if the converted file exits
                                        if not self.convertedCifToEMDSpaceExits(conv_f):
                                            # There is no converted file, convert it now
                                            convert(f, conv_f).em2emd()
                                    if os.path.exists(conv_f) and os.path.isfile(conv_f):
                                        of = os.path.join(self.outputXmlDir, id + '.xml')
                                        schema_loc = os.path.join(script_loc, self.schema)
                                        if os.path.exists(schema_loc):
                                            #print "conv_f %s" % conv_f
                                            #print "of %s" % of
                                            #print "schema_loc %s" % schema_loc
                                            translator.translate_and_validate(conv_f, of, schema_loc)
                                            #print "translation done"
                                            a_log = translator.current_entry_log
                                            log_id = a_log.id
                                        else:
                                            print 'schema missing. Add it here'
                                    else:
                                        print 'The file ' + conv_f + 'cannot be converted into the _emd space and therefore, cannot be translated'
                                else:
                                    # print "Nope"
                                    pass
                            except IOError as exp:
                                print exp
                            #print '\nPDBE LOGGER LOGS\n'
                            translator.write_logger_logs(write_error_log=True)
                #print '\nPDBE CLASS LOGS\n'
                if translator.is_translation_log_empty:
                    print 'NO ERRORS FOUND'
                else:
                    for entry_log in translator.translation_log.logs:
                        print entry_log.id
                        if entry_log.is_error_log_empty:
                            print "\nNO ERRORS FOUND FOR THIS ENTRY\n"
                        else:
                            for err in entry_log.errors:
                                print '\nERRORS: %s\n' % err.log_text
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
        latestEMDBDump = self.db_dump
        with open(latestEMDBDump) as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                ids.append(row[1])
                d_ids.append(row[0])
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
