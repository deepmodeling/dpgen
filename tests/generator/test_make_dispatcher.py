import os,sys,sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import make_dispatcher

class TestDispatcher(unittest.TestCase):
    # def test_ssh_slurm(self):
    #     dis = make_dispatcher({
    #         'batch':    'slurm',
    #         'hostname': 'localhost',
    #         'username': 'wanghan',
    #         'port':     22,
    #         'work_path': '.',
    #     })
    #     self.assertEqual(dis.context.__name__, 'SSHContext')
    #     self.assertEqual(dis.batch.__name__, 'Slurm')

    def test_local_slurm(self):
        dis = make_dispatcher({
            'batch':    'slurm',
            'work_path': '.',
        })
        self.assertEqual(dis.context.__name__, 'LocalContext')
        self.assertEqual(dis.batch.__name__, 'Slurm')

    def test_lazy_local_slurm(self):
        dis = make_dispatcher({
            'batch':    'slurm',
            'lazy_local': True,
            'work_path': '.',
        })
        self.assertEqual(dis.context.__name__, 'LazyLocalContext')
        self.assertEqual(dis.batch.__name__, 'Slurm')

    def test_dep_lazy_local_slurm(self):
        dis = make_dispatcher({
            'machine_type':    'slurm',
            'lazy_local': True,
            'work_path': '.',
        })
        self.assertEqual(dis.context.__name__, 'LazyLocalContext')
        self.assertEqual(dis.batch.__name__, 'Slurm')
