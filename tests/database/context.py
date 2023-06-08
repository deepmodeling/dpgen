import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))


def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
