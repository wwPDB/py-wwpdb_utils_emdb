__date__ = "2020-11-17"
__version__ = "0.1"


import os
from optparse import OptionParser

from wwpdb.utils.emdb.cif_emdb_translator.cif_emdb_translator import CifEMDBTranslator


def main():
    """
    This script calls CifEMDBTranslator to translate
    an mmcif emd file into an EMDB v3.x XML header file
    :return:
    """
    in_cif = None
    out_xml = None
    out_dir = "data/test/xml_v3_out/"
    schema = "emdb.xsd"
    usage = """
            test_cif_to_xml_translator.py [options]
            Convert an mmcif emd file to EMDB v3.x xml file

            Examples:
            python test_cif_to_xml_translator.py -i input_cif_file -o output_xml_file
            python test_cif_to_xml_translator.py -i "data/cif/EMD-0000.txt" -o "data/test/xml_v3_out/EMD-0000.xml"

            Typical run:
            python test_cif_to_xml_translator.py -i in.cif -o out.xml
                in.xml: a cif file in emd space
                out.cif: an EMDB XML file 3.x
            """
    parser = OptionParser(usage=usage, version=__version__)
    parser.add_option("-i", "--in-file", action="store", type="string", metavar="FILE", dest="in_cif", help="Input cif file")
    parser.add_option("-o", "--out-file", action="store", type="string", metavar="FILE", dest="out_xml", help="Output xml file")
    (options, args) = parser.parse_args()

    # if len(args) == 0:
    #     print(usage)
    #     sys.exit("No input options given!")

    if options.in_cif:
        in_cif = options.in_cif
    else:
        in_cif = "data/cif/EMD-0000.txt"
    if options.out_xml:
        out_xml = options.out_xml
    else:
        out_xml = out_dir + "EMD-0000.xml"

    # create the translator object
    translator = CifEMDBTranslator()
    # set up the translator's logger object
    translator.set_logger_logging(log_error=True, error_log_file_name=os.path.join(os.path.dirname(os.path.realpath(__file__)), "ERROR.log"))
    translator.set_logger_logging(True, True, True, False)

    # do the translation
    if os.path.exists(in_cif):
        if os.path.exists(out_dir):
            translator.translate_and_validate(in_cif, out_xml, schema)
            if os.path.exists(out_dir):
                print("Output XML file", out_xml, "has been created")
            translator.write_logger_logs(True, True, True)
        else:
            print("Output file path", out_xml, "does not exit")
    else:
        print("Input file", in_cif, "does not exist")


if __name__ == "__main__":
    main()
