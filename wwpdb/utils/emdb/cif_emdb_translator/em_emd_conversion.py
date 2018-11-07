import sys
import traceback, argparse

from wwpdb.api.plugins.UtilsBase import UtilsBase
from wwpdb.api.facade.ConfigInfo import ConfigInfo
from pdbx_v2.trans.InstanceMapper import InstanceMapper

class convert(UtilsBase):

    def __init__(self, inFile=None, outFile=None, verbose=False, log=sys.stderr):


        headerFilters = ['all', 'prereleasetitle']
        cI = ConfigInfo()
        self.mappingInfoPath = cI.get("SITE_EXT_DICT_MAP_EMD_FILE_PATH")
        self.inFile = inFile
        self.outFile = outFile
        self.verbose = verbose
        self.log = log
        self.im = InstanceMapper(verbose=self.verbose, log=self.log)
        self.im.setMappingFilePath(self.mappingInfoPath)
        #self.im.setFilterTagList(headerFilters)


    def em2emd(self):
        ok = self.im.translate(self.inFile, self.outFile, mode="src-dst")
        print ok

    def emd2em(self):
        ok = self.im.translate(self.inFile, self.outFile, mode="dst-src")
        print ok

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='em translate')
    parser.add_argument('-o', '--output', help='output file', type=str, required=True)
    parser.add_argument('-i', '--input', help='input file', type=str, required=True)
    parser.add_argument('-d', '--emToEmd', help='em to emd', action='store_true')
    parser.add_argument('-e', '--emdToEm', help='emd to em', action='store_true')

    args = parser.parse_args()

    if not sys.argv[1:]:
        parser.print_help()
        exit()

    if args.emToEmd:
        convert(inFile=args.input, outFile=args.output).em2emd()

    elif args.emdToEm:
        convert(inFile=args.input, outFile=args.output).emd2em()

    else:
        parser.print_help()
        exit()
