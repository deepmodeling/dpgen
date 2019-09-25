import sys,os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))

from dpgen.dispatcher.SSHContext import SSHSession
from dpgen.dispatcher.SSHContext import SSHContext
from dpgen.dispatcher.LocalContext import LocalSession
from dpgen.dispatcher.LocalContext import LocalContext
from dpgen.dispatcher.Shell import Shell
from dpgen.dispatcher.JobStatus import JobStatus
from dpgen.dispatcher.Dispatcher import Dispatcher

def my_file_cmp(test, f0, f1):
    with open(f0) as fp0 :
        with open(f1) as fp1:
            test.assertTrue(fp0.read() == fp1.read())

def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))

