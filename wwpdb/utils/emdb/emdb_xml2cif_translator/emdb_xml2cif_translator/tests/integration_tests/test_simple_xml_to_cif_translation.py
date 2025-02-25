import os
import emdb_xml2cif_translator.test_data.xml
import emdb_xml2cif_translator.test_data.cif
import emdb_xml2cif_translator.input_files
from emdb_xml2cif_translator.translator_classes.EMDBXmlToCifTranslator import EMDBXmlToCifTranslator
from emdb_xml2cif_translator.translator_classes.LoggingUtil import LoggingUtil


def main():
   
    translator = EMDBXmlToCifTranslator()
    translator.translate_xml_to_cif(os.path.join(emdb_xml2cif_translator.test_data.xml.__path__[0], "emd-0001.xml"),
                                    os.path.join(emdb_xml2cif_translator.test_data.cif.__path__[0], "emd-0001.cif"))


if __name__ == "__main__":
    usage = """
            Converts one EMDB header (xml) file (v3.x) to an mmcif file
            
            Example: python test_xml_to_cif_translator.py
            
            """

    main()
