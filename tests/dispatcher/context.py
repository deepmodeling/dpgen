import sys,os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from dpgen.dispatcher.LocalContext import LocalSession
from dpgen.dispatcher.LocalContext import LocalContext
from dpgen.dispatcher.LazyLocalContext import LazyLocalContext
from dpgen.dispatcher.SSHContext import SSHSession
from dpgen.dispatcher.SSHContext import SSHContext
# from dpgen.dispatcher.Dispatcher import FinRecord
from dpgen.dispatcher.Dispatcher import _split_tasks

from dpgen.dispatcher.LocalContext import _identical_files

def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
