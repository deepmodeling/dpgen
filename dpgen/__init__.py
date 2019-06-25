from   __future__ import unicode_literals, print_function
import logging
import os

# from . import lib

NAME="dpgen"
SHORT_CMD="dpgen"
dlog = logging.getLogger(__name__)
dlog.setLevel(logging.DEBUG)
dlogf = logging.FileHandler(os.getcwd()+os.sep+SHORT_CMD+'.log')
dlogf_formatter=logging.Formatter('%(asctime)s - %(levelname)s : %(message)s')
#dlogf_formatter=logging.Formatter('%(asctime)s - %(name)s - [%(filename)s:%(funcName)s - %(lineno)d ] - %(levelname)s \n %(message)s')
dlogf.setFormatter(dlogf_formatter)
dlog.addHandler(dlogf)

__author__    = "Han Wang"
__copyright__ = "Copyright 2019"
__version__   = "0.1.0"
__status__    = "Development"
__date__      = "Jun 27, 2019"




def info():
    """
        Show basic information about """+NAME+""", its location and version.
    """

    print('DeepModeling\n------------\n')
    print('Version: ' + __version__)
    print('Path:    ' + __path__[0])
    print('Date:    ' + __date__)
    print()

