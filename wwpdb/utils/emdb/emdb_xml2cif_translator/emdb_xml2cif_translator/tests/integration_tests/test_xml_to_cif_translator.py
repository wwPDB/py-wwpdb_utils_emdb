import os, argparse
from pathlib import Path
import emdb_xml2cif_translator.test_data.xml
import emdb_xml2cif_translator.test_data.cif
import emdb_xml2cif_translator.input_files
from emdb_xml2cif_translator.translator_classes.EMDBXmlToCifTranslator import EMDBXmlToCifTranslator
from emdb_xml2cif_translator.translator_classes.LoggingUtil import LoggingUtil

def main():
    xml_input_dir = os.path.join(emdb_xml2cif_translator.test_data.xml.__path__[0])
    autodep_input_dir = os.path.join(emdb_xml2cif_translator.test_data.xml.__path__[0])
    cif_output_dir = os.path.join(emdb_xml2cif_translator.test_data.cif.__path__[0])

    # Specify the XML and autodep input filenames
    xml_files = [f for f in os.listdir(xml_input_dir) if f.startswith('emd-')]
    autodep_files = [f for f in os.listdir(xml_input_dir) if f.startswith('autodep')]

    # Create translation_list with specified input files
    translation_list = {
        os.path.join(xml_input_dir, xml_file): [
            os.path.join(autodep_input_dir, autodep_file),  # Second input file (autodep_input_dir)
            os.path.join(cif_output_dir, f"{os.path.splitext(xml_file)[0].rsplit('-', 1)[0]}.cif")
        ] for xml_file, autodep_file in zip(xml_files, autodep_files)
    }

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

    mapping_file = os.path.join(emdb_xml2cif_translator.input_files.__path__[0], "sample-mappings.txt")

    for xml_input_file, (autodep_input_file, cif_output_file) in translation_list.items():
        translator = EMDBXmlToCifTranslator(logging_utility)
        # translator = EMDBXmlToCifTranslator(logging_utility, mapping_file)
        translator.translate_xml_to_cif(xml_input_file, autodep_input_file, cif_output_file)
        # translator.output_translation_log()

if __name__ == "__main__":
    prog = "XML2CIF"
    usage = """
            Converting legacy entries XML files to internal mmcif format
            Example:
            python test_xml_to_cif_translator.py 
          """
    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    # parser.add_argument('-x', '--xmlDir', type=Path, help="Directory path to the legacy entries XML header files .")
    # parser.add_argument('-c', '--cifDir', type=Path, help="Directory path to write the output mmcif files.")
    args = parser.parse_args()

    main()

