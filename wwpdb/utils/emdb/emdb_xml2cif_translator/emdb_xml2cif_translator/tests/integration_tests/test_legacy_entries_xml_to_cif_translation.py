import os, argparse
from glob import glob
import emdb_xml2cif_translator.test_data.xml
import emdb_xml2cif_translator.test_data.cif
import emdb_xml2cif_translator.input_files
from emdb_xml2cif_translator.translator_classes.EMDBXmlToCifTranslator import EMDBXmlToCifTranslator
from emdb_xml2cif_translator.translator_classes.LoggingUtil import LoggingUtil
from emdb_xml2cif_translator.translator_classes.LegacyEntries import LegacyEntries

def main():
    legacy_set = LegacyEntries().ids
    
    base_directory = os.path.join(emdb_xml2cif_translator.test_data.xml.__path__[0])
    legacy_folder = glob(os.path.join(base_directory, 'EMD-*'))
    for file in legacy_folder:
        id_num = file.rsplit('/', 1)[1].split('-', 1)[1]
        print(f"Translating EMD-{id_num}")
        xml_filepath = os.path.join(file, f"structures/EMD-{id_num}/header/emd-{id_num}.xml")
        autodep_filepath = os.path.join(file, f"autodep.xml")
        cif_output_dir = os.path.join(emdb_xml2cif_translator.test_data.cif.__path__[0])
        if not os.path.isfile(xml_filepath):
            print(f"{xml_filepath} not found.")
            return None
        # Create translation_list with specified input files
        translation_list = {
            xml_filepath: [
                autodep_filepath,
                os.path.join(cif_output_dir, f"emd-{id_num}.cif")
            ]}

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
            Converting entries XML files to mmcif format
            By default, the output will be in internal mmcif format.
            Example:
            python test_xml_to_cif_translator.py 
          """
    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-p', '--public_mmcif_format', type=bool, help="Flag to indicate if the output should be in public mmcif format.")
    parser.add_argument('-l', '--legacy', type=bool, default=False, help="Flag to indicate if the translation is .")
    # parser.add_argument('-x', '--xmlDir', type=Path, help="Directory path to the legacy entries XML header files .")
    # parser.add_argument('-c', '--cifDir', type=Path, help="Directory path to write the output mmcif files.")
    args = parser.parse_args()

    main()

