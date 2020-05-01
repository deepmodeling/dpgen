import os,sys
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import parse_cur_job
from .context import param_file
from .context import machine_file
from .context import setUpModule

class TestParseCurJob(unittest.TestCase):
    def test_npt (self) :        
        ens = 'npt'
        ts = [100,200]
        ps = [1e5,1e6,1e7]
        ns = 1000
        tf = 10
        cur_job = {}
        cur_job['ens'] = ens
        cur_job['Ts'] = ts
        cur_job['Ps'] = ps
        cur_job['nsteps'] = ns
        cur_job['t_freq'] = tf
        res = parse_cur_job(cur_job)
        for ii,jj in zip(res, [ens, ns, tf, ts, ps, None, None]) :
            self.assertEqual(ii,jj)

    def test_nvt (self) :        
        ens = 'nvt'
        ts = [100,200]
        ps = [1e5,1e6,1e7]
        ns = 1000
        tf = 10
        cur_job = {}
        cur_job['ens'] = ens
        cur_job['Ts'] = ts
        cur_job['Ps'] = ps
        cur_job['nsteps'] = ns
        cur_job['t_freq'] = tf
        res = parse_cur_job(cur_job)
        for ii,jj in zip(res, [ens, ns, tf, ts, [-1], None, None]) :
            self.assertEqual(ii,jj)

    def test_pka (self) :        
        ens = 'nvt'
        ts = [100,200]
        ps = [1e5,1e6,1e7]
        ns = 1000
        tf = 10
        pka = [10, 20, 30]
        dt = 0.001
        cur_job = {}
        cur_job['ens'] = ens
        cur_job['Ts'] = ts
        cur_job['Ps'] = ps
        cur_job['nsteps'] = ns
        cur_job['t_freq'] = tf
        cur_job['pka_e'] = pka
        cur_job['dt'] = dt
        res = parse_cur_job(cur_job)
        for ii,jj in zip(res, [ens, ns, tf, ts, [-1], pka, dt]) :
            self.assertEqual(ii,jj)


if __name__ == '__main__':
    unittest.main()

