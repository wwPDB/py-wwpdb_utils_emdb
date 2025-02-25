from re import I
import sys
from emdb_xml2cif_translator.translator_classes.CIF import CIF
from emdb_xml2cif_translator.translator_classes.EMDBMetadata import EMDBMetadata
from emdb_xml2cif_translator.translator_classes.LoggingUtil import LoggingUtil
import os
import errno

FILE_NF_ERROR = Exception
if sys.version_info.major == 3:
    FILE_NF_ERROR = FileNotFoundError


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
        self.__autodep_file_read = False
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
                raise FILE_NF_ERROR(errno.ENOENT, os.strerror(errno.ENOENT), input_file)
            else:
                self.data.emdb_id_from_filename = os.path.splitext(os.path.basename(input_file))[0]
                self.__emdb_header_file_read = self.data.parse_into_xml_tree(input_file)

    def read_autodep_file(self, autodep_input_file=None):
        """
        This method reads the autodep input file and saves its values into the EMDBMetaData object as a dictionary
        If autodep_input_file is either not given or doesn't exist, an exception is raised.
        The autodep input file data is stored into the EMDBMetadata object (self.data)

        :param autodep_input_file: An XML file; autodep input file
        :return __autodep_file_read: True if the XML is successfully parsed
        """
        autodep_input_file = self.get_autodep_file()
        if not autodep_input_file:
            raise Exception("Translator cannot read the autodep.xml file as the file name is not given")
        else:
            if not os.path.exists(autodep_input_file):
                raise FILE_NF_ERROR(errno.ENOENT, os.strerror(errno.ENOENT), autodep_input_file)
            else:
                self.__autodep_file_read = self.data.parse_autodep(autodep_input_file)

    def get_autodep_file(self):
        """
        This method returns the autodep input file name

        :return autodep_input_file: a string; autodep input file name
        """
        autodep_input_file = None
        if self.data:
            autodep_input_file = self.data.get_autodep_file()
        return autodep_input_file
    
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

    def read_input(self, input_xml_file):
        """
        The method reads the input XML file and the autodep file if the input XML file is a legacy entry

        :param input_xml_file (_type_): an XML file; EMDB entry header file; a v3.x EMDB file
        """
        if not self.__emdb_header_file_read and self.data:
            # read input file now and save the result in self.data.emd
            self.read_emdb_header_file(input_xml_file)
            if self.data.is_legacy:
                # read autodep input file and save the result in self.data.autodep
                self.read_autodep_file()

    def set_logging(self):
        """
        The method sets the logging framework for the translation process
        """
        # set the logging framework
        if self.__emdb_header_file_read and not self.logging:
            # set the logging lists
            if self.data.emdb_id_from_filename:
                self.translation_log = LoggingUtil.EntryLogs(self.data.emdb_id_from_filename)
                if self.logging and self.translation_log:
                    self.__logger_set = True

    def set_output(self, output_cif_file):
        """
        The method sets the output cif file name and creates a CIF object within the EMDBData object

        :param output_cif_file: a cif file
        """
        if not self.__emdb_header_file_translated and not self.__autodep_file_translated:
            if not output_cif_file:
                # If the output file is not give a default output file name is used;
                # output file is always written
                output_cif_file = self.DEFAULT_OUTPUT_FILENAME
            # Create a CIF object within the EMDBData object to write the output mmcif file
            self.data.cif = CIF(output_cif_file)

    def translate(self):
        """
        The method calls the process method of the EMDBMetadata object to translate the input data
        """
        if self.__emdb_header_file_read and self.data:
            # translate the data now
            return self.data.process() 

    def translate_xml_to_cif(self, input_xml_file, autodep_xml_file=None, output_cif_file=None):
        """
        The method reads the input file(s) and writes the result of a xml-to-cif translation of data from
        'input_xml_file' into 'output_mmcif_file'

        :param input_xml_file: an XML file; EMDB entry header file; a v3.x EMDB file
        :param output_cif_file: a cif file
        :return translated_successfully: a boolean; True is all translation components are performed successfully
        """
        # before translating starts, the input header and if a legacy entry the autodep file need to be read;
        # at this point the files are still not read and the read flags should not be set, however
        # the data object has been created and initialised when the translator object was created
        self.read_input(input_xml_file, autodep_xml_file)
        # input is read, EMDB metadata object created. The logger needs to be established 
        self.set_logging()       
        # The logger is established. 
        # The output cif file name and the cif object needs creating
        self.set_output(output_cif_file)
        # All set - translate now
        if self.translate():
            # translation is done, write the output cif file
            self.write_mmcif_file()