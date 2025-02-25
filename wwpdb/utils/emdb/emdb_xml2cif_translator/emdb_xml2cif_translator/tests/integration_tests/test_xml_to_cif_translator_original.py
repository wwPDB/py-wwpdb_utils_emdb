import os
import emdb_xml2cif_translator.test_data.xml
import emdb_xml2cif_translator.test_data.cif
import emdb_xml2cif_translator.input_files
from emdb_xml2cif_translator.translator_classes.EMDBXmlToCifTranslator import EMDBXmlToCifTranslator
from emdb_xml2cif_translator.translator_classes.LoggingUtil import LoggingUtil


def main():
   
    xml_input_file_0001 = os.path.join(emdb_xml2cif_translator.test_data.xml.__path__[0], "emd-0001.xml")

    cif_output_file_0001 = os.path.join(emdb_xml2cif_translator.test_data.cif.__path__[0], "emd-0001.cif")

    translation_list = {xml_input_file_0001: cif_output_file_0001}

    logging_params = {
        "info": {
            "log_file": {
                "log": True,
                "name": "info.log"
            },
            "log_stream": {
                "log": True,
                "name": None
            }
        },
        "warn": {
            "log_file": {
                "log": True,
                "name": "warn.log"
            },
            "log_stream": {
                "log": True,
                "name": None
            }
        },
        "error": {
            "log_file": {
                "log": True,
                "name": "error.log"
            },
            "log_stream": {
                "log": True,
                "name": None
            }
        },
        "console": True
    }

    logging_utility = LoggingUtil(logging_params) 

    for xml_input_file, cif_output_file in translation_list.items():
        translator = EMDBXmlToCifTranslator(logging_utility)
        translator.translate_xml_to_cif(xml_input_file, cif_output_file)
        # translator.output_translation_log()


if __name__ == "__main__":
    usage = """
            test_xml_to_cif_translator.py [options]
            
            Converts one EMDB header (xml) file (v3.x) to an mmcif file
            
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

    main()
