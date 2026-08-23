import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "generator"
from dpgen.generator.arginfo import model_devi_jobs_args
from dpgen.generator.run import _get_lammps_job_settings

from .context import (
    parse_cur_job,
    setUpModule,  # noqa: F401
)


class TestParseCurJob(unittest.TestCase):
    def test_npt(self):
        ens = "npt"
        ts = [100, 200]
        ps = [1e5, 1e6, 1e7]
        ns = 1000
        tf = 10
        cur_job = {}
        cur_job["ens"] = ens
        cur_job["Ts"] = ts
        cur_job["Ps"] = ps
        cur_job["nsteps"] = ns
        cur_job["t_freq"] = tf
        res = parse_cur_job(cur_job)
        for ii, jj in zip(res, [ens, ns, tf, ts, ps, None, None]):
            self.assertEqual(ii, jj)

    def test_nvt(self):
        ens = "nvt"
        ts = [100, 200]
        ps = [1e5, 1e6, 1e7]
        ns = 1000
        tf = 10
        cur_job = {}
        cur_job["ens"] = ens
        cur_job["Ts"] = ts
        cur_job["Ps"] = ps
        cur_job["nsteps"] = ns
        cur_job["t_freq"] = tf
        res = parse_cur_job(cur_job)
        for ii, jj in zip(res, [ens, ns, tf, ts, [-1], None, None]):
            self.assertEqual(ii, jj)

    def test_pka(self):
        ens = "nvt"
        ts = [100, 200]
        ps = [1e5, 1e6, 1e7]
        ns = 1000
        tf = 10
        pka = [10, 20, 30]
        dt = 0.001
        cur_job = {}
        cur_job["ens"] = ens
        cur_job["Ts"] = ts
        cur_job["Ps"] = ps
        cur_job["nsteps"] = ns
        cur_job["t_freq"] = tf
        cur_job["pka_e"] = pka
        cur_job["dt"] = dt
        res = parse_cur_job(cur_job)
        for ii, jj in zip(res, [ens, ns, tf, ts, [-1], pka, dt]):
            self.assertEqual(ii, jj)

    def test_job_local_lammps_settings_override_global_defaults(self):
        cur_job = {"dt": 0.001, "neidelay": 5, "taut": 0.2, "taup": 1.0}
        jdata = {
            "model_devi_dt": 0.002,
            "model_devi_neidelay": 10,
            "model_devi_taut": 0.1,
            "model_devi_taup": 0.5,
        }
        self.assertEqual(_get_lammps_job_settings(cur_job, jdata), (0.001, 5, 0.2, 1.0))

    def test_job_schema_accepts_local_timestep(self):
        arginfo = model_devi_jobs_args()
        jobs = [
            {
                "sys_idx": [0],
                "ensemble": "nvt",
                "temps": [300.0],
                "nsteps": 100,
                "trj_freq": 10,
                "dt": 0.001,
            }
        ]
        normalized = arginfo.normalize_value(jobs, trim_pattern="_*")
        arginfo.check_value(normalized, strict=True)


if __name__ == "__main__":
    unittest.main()
