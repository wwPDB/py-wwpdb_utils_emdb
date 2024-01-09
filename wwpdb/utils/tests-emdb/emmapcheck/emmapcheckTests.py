import os
import django

os.environ["DJANGO_SETTINGS_MODULE"] = "wwpdb.apps.deposit.settings"
django.setup()

from wwpdb.apps.deposit.depui.FileConversion import FileConversion
from wwpdb.apps.deposit.depui.common_functions import getFile
from wwpdb.utils.wf.dbapi.WfDbApi import WfDbApi

username = "D_800288"

wf_op = "emmapcheck_dep"
fc = FileConversion()
wfApi = WfDbApi(verbose=True)
status = fc.wfStart(username, wfApi, wf_op, "wf_op_emmapcheck_fs_tempdep.xml")
wf_data = fc.wfStatus(username, wfApi)
finished = fc.wfFinished(username, wfApi, wf_data, wf_op)

if status == 0 or not finished:
    print("Unable to perform checks on uploaded files!")

path2output = getFile(username, "tempdep", "em-map-report", "json", versionId='latest', mileStone=None, wfInstanceId=None, partNumber=None)


# import sys
# if __name__ == '__main__':
#     # unittest.main()
#     sys.exit(0)