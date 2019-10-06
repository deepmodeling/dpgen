from   __future__ import unicode_literals, print_function
import logging
import os


ROOT_PATH=__path__[0]
NAME="dpgen"
SHORT_CMD="dpgen"
dlog = logging.getLogger(__name__)
dlog.setLevel(logging.INFO)
dlogf = logging.FileHandler(os.getcwd()+os.sep+SHORT_CMD+'.log')
dlogf_formatter=logging.Formatter('%(asctime)s - %(levelname)s : %(message)s')
#dlogf_formatter=logging.Formatter('%(asctime)s - %(name)s - [%(filename)s:%(funcName)s - %(lineno)d ] - %(levelname)s \n %(message)s')
dlogf.setFormatter(dlogf_formatter)
dlog.addHandler(dlogf)

__author__    = "Han Wang"
__copyright__ = "Copyright 2019"
__version__   = "0.4.0"
__status__    = "Development"
__date__      = "Nov 5, 2019"



def info():
    """
        Show basic information about """+NAME+""", its location and version.
    """

    print('DeepModeling\n------------')
    print('Version: ' + __version__)
    print('Date:    ' + __date__)
    print('Path:    ' + ROOT_PATH)
    print('')
    print('Dependency')
    print('------------')
    for modui in ['numpy', 'dpdata', 'pymatgen', 'monty', 'ase', 'paramiko', 'custodian' ]:
        try:
            mm = __import__(modui)
            print('%10s %10s   %s' % (modui, mm.__version__, mm.__path__[0]))
        except ImportError:
            print('%10s %10s Not Found' % (modui, ''))
    print()
