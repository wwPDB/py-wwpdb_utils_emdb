import sys
from emdb_xml2cif_translator.translator_classes.CIF import CIF
from emdb_xml2cif_translator.translator_classes.EMDBMetadata import EMDBMetadata
from emdb_xml2cif_translator.translator_classes.LoggingUtil import LoggingUtil
import os
import errno

file_nf_error = Exception
if sys.version_info.major == 3:
    file_nf_error = FileNotFoundError


class EMDBXmlToCifTranslator(object):
    """
    This class provides functionality for translating a given EMDB v3.x XML header file into a cif output file
    When initialising, the class object can accept an LoggingUtil object for the translation logging purposes
    """
    # default output cif file name in case one is not provided
    DEFAULT_OUTPUT_FILENAME = 'output.cif'

    def __init__(self, logging_util=None):
        # flags for indicating a translation in process
        self.__emdb_header_file_read = False
        self.__emdb_header_file_translated = False
        self.__mmcif_file_written = False
        self.__mmcif_file_validates = False
        self.__logger_set = False

        # objects aiding a translation
        # logging is LoggingUtil object
        self.logging = logging_util
        # translation_log is an EntryLogs object
        self.translation_log = None
        # data is an object of the EMDBMetaData class containing:
        # 1. data.map_utils: a dictionary containing xml-to-cif mappings logic
        # 2. data.emd: a dictionary containing values from the input XML file
        # 3. data.cif_ready_data: a dictionary for data.cif, an object of the CIF() class
        self.data = EMDBMetadata()

    def read_emdb_header_file(self, input_file):
        """
        This method reads the input file and saves its values into the EMDBMetaData object as a dictionary
        If input_file is either not given or doesn't exist an exception is raised
        The input file data is stored into the EMDBMetadata object (self.data)

        :param input_file: An XML file; EMDB entry v3.x header file
        :return __emdb_header_file_read: True if the XML is successfully parsed
        """
        if not input_file:
            raise Exception("Translator cannot read the input file as the file name is not given")
        else:
            if not os.path.exists(input_file):
                raise file_nf_error(errno.ENOENT, os.strerror(errno.ENOENT), input_file)
            else:
                self.data.emdb_id_from_filename = os.path.splitext(os.path.basename(input_file))[0]
                self.__emdb_header_file_read = self.data.parse_into_xml_tree(input_file)

        return self.__emdb_header_file_read

    def write_mmcif_file(self):
        """
        This method checks if the translation is performed and calls onto the cif object to write the output file

        :return __mmcif_file_written: a boolean; True if the cif object has successfully written cif file
        """
        if self.__emdb_header_file_translated:
            if not self.data.cif.filename:
                raise Exception("Translator cannot write output file as the file name is not given")
            elif not self.__mmcif_file_written:
                self.__mmcif_file_written = self.data.cif.write()

        return self.__mmcif_file_written

    def translate_xml_to_cif(self, input_xml_file, output_cif_file=None):
        """
        The method reads 'input_xml_file' and writes the result of a xml-to-cif translation of data from
        'input_xml_file' into 'output_mmcif_file'

        :param input_xml_file: an XML file; EMDB entry header file; a v3.x EMDB file
        :param output_cif_file: a cif file
        :return translated_successfully: a boolean; True is all translation components are performed successfully
        """
        translated_successfully = False

        # before translating starts, the input file needs to be read;
        # at this point the input file is still not read and the flag should not be set
        if not self.__emdb_header_file_read:
            # data object has been created and initialised when the translator object was created
            if self.data:
                # read input file now and save the result in self.data.emd
                if self.read_emdb_header_file(input_xml_file):
                    # set the logging framework
                    if not self.logging:
                        # set the logging lists
                        if self.data.emdb_id_from_filename:
                            self.translation_log = LoggingUtil.EntryLogs(self.data.emdb_id_from_filename)
                            if self.logging and self.translation_log:
                                self.__logger_set = True
                    if not self.__emdb_header_file_translated:
                        # Already set for a translation: input XML file is read, EMDB metadata object created and
                        # the logger established
                        # Still needed: output filename and cif object
                        if not output_cif_file:
                            # If the output file is not give a default output file name is used;
                            # output file is always written
                            output_cif_file = self.DEFAULT_OUTPUT_FILENAME
                        # Create a CIF object within the EMDBData object to write the output mmcif file
                        self.data.cif = CIF(output_cif_file)
                        # All set - translate now
                        self.__emdb_header_file_translated = self.data.process()
                    if self.__emdb_header_file_translated:
                        if self.write_mmcif_file():
                            translated_successfully = True

        return translated_successfully
