
import os
import test_data.xml
import test_data.cif
import input_files
from translator_classes.EMDBXmlToCifTranslator import EMDBXmlToCifTranslator
from translator_classes.LoggingUtil import LoggingUtil


def main():
    """
    This script calls EMDBXmlToCifTranslator to translate
    an EMDB v3.x XML header file into an mmcif file
    """
    usage = """
            test_xml_to_cif_translator.py [options]
            Convert a EMDB v3.x xml file to mmcif file
            
            Examples: 
            python test_xml_to_cif_translator.py -i input_xml_file -o output_cif_file
            
            Typical run:
            python test_xml_to_cif_translator.py -i in.xml -o out.cif
                in.xml: an EMDB XML file 3.x
                out.cif: a cif file in emd space
            
            python test_xml_to_cif_translator.py -i in.xml -o out.cif -d in.dic
                 in.dic: An mmcif dictionary that the out.cif is translated from
            """
    # print(usage)

    xml_input_file_0001 = os.path.join(test_data.xml.__path__[0], "emd-0001.xml")

    cif_output_file_0001 = os.path.join(test_data.cif.__path__[0], "emd-0001.cif")

    translation_list = {xml_input_file_0001: cif_output_file_0001}

    logging_params = {
        "info": {
            "log_file": {
                "log": True,
                "name": "/Users/sanja/IdeaProjects/emdb-xml2cif-translator/logs/info.log"
            },
            "log_stream": {
                "log": True,
                "name": None
            }
        },
        "warn": {
            "log_file": {
                "log": True,
                "name": "/Users/sanja/IdeaProjects/emdb-xml2cif-translator/logs/warn.log"
            },
            "log_stream": {
                "log": True,
                "name": None
            }
        },
        "error": {
            "log_file": {
                "log": True,
                "name": "/Users/sanja/IdeaProjects/emdb-xml2cif-translator/logs/error.log"
            },
            "log_stream": {
                "log": True,
                "name": None
            }
        },
        "console": True
    }

    logging_utility = LoggingUtil(logging_params)

    mapping_file = os.path.join(input_files.__path__[0], "sample-mappings.txt")

    for xml_input_file, cif_output_file in translation_list.items():
        translator = EMDBXmlToCifTranslator(logging_utility)
        # translator = EMDBXmlToCifTranslator(logging_utility, mapping_file)
        translator.translate_xml_to_cif(xml_input_file, cif_output_file)
        # translator.output_translation_log()


if __name__ == "__main__":
    main()
