import os
from dpgen import dlog
from pymatgen.ext.matproj import MPRester, MPRestError
from pymatgen.core import Structure

web="materials.org"

def check_apikey():
    try:
      apikey=os.environ['MAPI_KEY']
    except KeyError:
       print("You have to get a MAPI_KEY from "+web)
       print("and execute following command:")
       print('echo "export MAPI_KEY=yourkey">> ~/.bashrc')
       print("source ~/.bashrc")
       os._exit(0)
    try:
      return MPRester(apikey)
    except MPRestError:
       dlog.info("MPRester Error, you need to prepare POSCAR manually")
       os._exit(0)

def get_structure(mp_id):
    mpr=check_apikey()
    return mpr.get_structure_by_material_id(mp_id)
