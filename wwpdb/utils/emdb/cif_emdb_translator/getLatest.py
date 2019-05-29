import re, glob, os, argparse, subprocess, shlex
from sqlalchemy import create_engine


SOFTWARE_ROOT= '/ebi/msd/work2/pdbdep/software'
PYMOL_CMD= os.path.join(SOFTWARE_ROOT,'bin', 'pymol_exe')
SVN_PATH= os.path.join(SOFTWARE_ROOT, 'scripts/image_scripts/production/2015onwards/trunk/images')
SCRIPT_PATH = os.path.join(SVN_PATH, 'PDB_plugin.py')
PREPARE= '/nfs/msd/oldnas9a/msd-pdbprepare/Processing/prepare'
EM_PREPARE = '/nfs/msd/oldnas9a/msd-emprepare/prepare/'
PRE_FTP = '/nfs/msd/work2/ftp/pdb/data/structures/all/mmCIF/'

getLatestPath = os.path.join(SOFTWARE_ROOT, 'scripts/cif_utilities/getPDB/getLatest.py')


# Database connection details
db_user = 'pdberead'
db_pass = 'pdberead55'
db_sid  = 'mysql-pdbe-onedep-prod.ebi.ac.uk'
db_port = "4437"

# Create a SQL Alchemy 'engine', then get the connection object so we can use normal SQL rather than objects
db_connection_string = "mysql://"+db_user+":"+db_pass+"@"+db_sid+":"+db_port

dataPath = '/nfs/msd/services/onedep/data/production/archive'

#d_and_a_machine = 'triton-2'
d_and_a_machine = 'pdbe-onedep-emdb'


class getWhat():
    def __init__(self, entry):
        self.entry = entry

    def getSingle(self):
        command = "ssh %s 'python %s %s'" %(d_and_a_machine, getLatestPath, self.entry)
        #print command
        args = shlex.split(command)
        proc = subprocess.Popen(args, stdout=subprocess.PIPE)
        fileType, location, copyStatic = re.split(',', proc.communicate()[0])
        if location:
            if copyStatic.strip() == 'copy':
                status, location = copyFile(location)
                if status != 0:
                    print 'copy failed'
                    return None, None
            return fileType, location
        else:
            return None, None

    def getPDB(self):
        command = "ssh %s 'python %s -p %s'" %(d_and_a_machine, getLatestPath, self.entry)
        #print command
        args = shlex.split(command)
        proc = subprocess.Popen(args, stdout=subprocess.PIPE)
        location = proc.communicate()[0]
        if location:
            status, location = copyFile(location)
            if status != 0:
                print 'copy failed'
                return None
            return location
        else:
            return None

    def getMapModel(self):
        command = "ssh %s 'python %s -e %s'" %(d_and_a_machine, getLatestPath, self.entry)
        #print command
        args = shlex.split(command)
        proc = subprocess.Popen(args, stdout=subprocess.PIPE)
        result = proc.communicate()[0]
        #print result
        if result:
            mapFile, cifFile, imageFile = re.split(',', result)

            if mapFile.strip() not in [None, 'None']:
                status, mapFile = copyFile(mapFile.strip())
                if status != 0:
                    print 'copy failed'
            if cifFile.strip() not in [None, 'None']:
                status, cifFile = copyFile(cifFile.strip())
                if status != 0:
                    print 'copy failed'

            if imageFile.strip() not in [None, 'None']:
                status, imageFile = copyFile(imageFile.strip())
                if status != 0:
                    print 'copy failed'

            return mapFile, cifFile, imageFile
        else:
            return None, None, None


def copyFile(filePath):

    if filePath:
        filePath = filePath.strip()
        user = os.getenv('USER')
        basename = os.path.basename(filePath)
        destination = os.path.join(os.getcwd(), basename)
        source = '%s@%s:%s' %(user,d_and_a_machine, filePath)
        command = 'scp %s %s' %(source, destination)
        print command
        #os.system(command)
        #return 0, basename

        p = subprocess.check_call(["scp", source, destination])
        if p == 0:
            return 0, basename
        else:
            return 1, 'copy failed'
    return 0, None

class PDBprocessedWhere():

    def __init__(self, file):
        self.file = file

    def Extension(self):
        shortName = os.path.basename(self.file)
        extension = os.path.splitext(shortName)[1]

        if extension:
            #print 'has extension of %s' %extension
            if extension in ['.pdb', '.ent']:
                return 'pdb', self.file, 'static'
            elif extension in ['.cif']:
                return 'cif', self.file, 'static'
            else:
                return 'cif', self.file, 'static'
        else:
            #print 'no extension'
            lower_entry = self.file.lower()
            ftp_cif = os.path.join(PRE_FTP, lower_entry + '.cif.gz')
            # check for entries which have been migrated to deppy
            entry = os.path.join('/ebi/msd/work2/w3_pdb05/MIGRATION', lower_entry)
            #print entry
            if os.path.exists(os.path.join('/ebi/msd/work2/w3_pdb05/MIGRATION', lower_entry)):
                #print 'yes'
                pdb_folder = os.path.join(PREPARE, lower_entry)
                ebi_cif = os.path.join(pdb_folder, 'ebi', lower_entry + '.cif')
                ebi_pdb = os.path.join(pdb_folder, 'ebi', 'pdb'+lower_entry+'.ent')
                orig_pdb = os.path.join(pdb_folder, 'ORIG', 'pdb'+lower_entry+'.ent')

                if os.path.exists(pdb_folder):
                    #print 'one of our entries'
                    if os.path.exists(ebi_cif):
                        return 'cif', ebi_cif, 'static'
                    elif os.path.exists(ebi_pdb):
                        return 'pdb', ebi_pdb, 'static'
                    elif os.path.exists(orig_pdb):
                        return 'pdb', orig_pdb, 'static'
                    else:
                        return None, None

                else:
                    depID = getDepID(lower_entry).checkType()
                    if depID:
                        # print depID
                        modelCif = getLatest(depID=depID, cifType='model', extension='cif').getFileList()
                        return 'cif', modelCif, 'copy'
                    elif os.path.exists(ftp_cif):
                        return 'cif', ftp_cif, 'static'
                    else:
                        return None, None
            else:
                depID = getDepID(lower_entry).checkType()
                #print " here"
                if depID:
                    #print depID
                    modelCif = getLatest(depID=depID, cifType='model', extension='cif').getFileList()
                    return 'cif', modelCif, 'copy'
                else:
                    return None, None


class getLatest():

    def __init__(self, depID, cifType, extension):

        self.depID = depID
        self.cifType = cifType
        self.extension = extension

    def getFileList(self):
        archivePath = os.path.join(dataPath, self.depID)
        #print archivePath
        filePattern = os.path.join(archivePath, self.depID+"_"+self.cifType+"_P1."+self.extension+".V")
        #print filePattern
        file_list = glob.glob(filePattern + "*")
        #print file_list
        if file_list:
            #print file_list
            fileName = self.prepare_file_list(file_list)
            return fileName
        else:
            return None

    def prepare_file_list(self, file_list):
        file_name = None
        if file_list:
            file_num = 0
            for f in file_list:
                if not '-cif-parser' in f:
                    num = int(re.split("\.", re.split(".V", f)[-1])[0])
                    if num > file_num:
                        file_num = num
                        file_name = f
        return file_name



class getDepID():

    def __init__(self, entryID):

        engine = create_engine(db_connection_string + '/status')
        self.connection = engine.connect()
        self.entryID = entryID
        self.depID = None
        self.pdbID = None
        self.emdbID = None

    def checkType(self):
        #print self.entryID
        if self.entryID[0:2] == "d_":
            self.depID = self.entryID.upper()
        elif self.entryID[0:3].lower() == 'emd':
            #print 'looks like emd'
            self.emQuery()
        else:
            self.sqlAccessionQuery()
        if self.depID:
            #print self.depID
            return self.depID
        else:
            #print 'no depID'
            return None

    def sqlAccessionQuery(self):

        self.pdbQuery()
        self.emQuery()

    def pdbQuery(self):
        #try to see if it is a PDB ID
        #print 'trying to see if %s is a pdbid' %self.entryID.upper()
        db_query = self.connection.execute("select * from dep_last_instance where pdb_id = '%s';" %self.entryID.upper()).fetchall()
        if db_query:
            for x in db_query:
                self.depID = x.dep_set_id
                self.pdbID = self.entryID
                self.emdbID = x.dep_emdb_id
        #print self.depID

    def emQuery(self):
        #print 'try to see if it is an EMDB ID'
        if self.entryID[0:3].lower() != 'emd':
            self.entryID = 'emd-%s' % self.entryID

        db_query = self.connection.execute("select * from dep_last_instance where dep_emdb_id = '%s';" % self.entryID.upper()).fetchall()
        if db_query:
            for x in db_query:
                self.depID = x.dep_set_id
                self.emdbID = self.entryID
                self.pdbID = x.pdb_id
        #print self.depID

class getMethod():

    def __init__(self, depID):
        engine = create_engine(db_connection_string + '/status')
        self.connection = engine.connect()
        self.depID = depID

    def getMethod(self):
        db_query = self.connection.execute("select * from dep_last_instance where dep_set_id = '%s';" % self.depID.upper()).fetchall()
        if db_query:
            for x in db_query:
                self.expMethodList = re.split(',' , x.dep_exp_method)
        return self.expMethodList

def returnEM(entry):
    depID = getDepID(entry).checkType()
    #print depID
    if depID:
        #check method type
        methodList = getMethod(depID=depID).getMethod()
        mapFile = getLatest(depID=depID, cifType='em-volume', extension='map').getFileList()
        cifFile = getLatest(depID=depID, cifType='model', extension='cif').getFileList()
        imageFile = getLatest(depID=depID, cifType='img-emdb', extension='png').getFileList()
        return mapFile, cifFile, imageFile
        #print '%s,%s' %(cifFile, mapFile)
    else:
        #print 'Not D&A'
        return None, None

def returnPDB(entry):
    depID = getDepID(entry).checkType()
    #print depID
    if depID:
        #check method type
        methodList = getMethod(depID=depID).getMethod()
        pdbFile = getLatest(depID=depID, cifType='model', extension='pdb').getFileList()

        return pdbFile
    else:
        return None


if '__main__' in __name__:
    parser = argparse.ArgumentParser(prog='getLatest.py')
    parser.add_argument('entry', help='entry to get latest file', type=str)
    parser.add_argument('-e', '--em', help='get em map and model', action="store_true")
    parser.add_argument('-p', '--pdb', help='get latest pdb', action="store_true")
    args = parser.parse_args()

    if args.em:
        mapFile, cifFile, imageFile = returnEM(args.entry)
        print '%s,%s,%s' %(mapFile, cifFile, imageFile)

    elif args.pdb:
        pdbFile = returnPDB(args.entry)
        print '%s' %(pdbFile)

    elif args.entry:
        fileType, location, copyStatic = PDBprocessedWhere(args.entry).Extension()
        print '%s,%s,%s' %(fileType, location, copyStatic)

        #depID = getDepID(args.entry).checkType()
        #modelCif = getLatest(depID=depID, cifType='model', extension='cif').getFileList()
        #print modelCif

    else:
        parser.print_help()
        exit()
