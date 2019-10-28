import os,dpdata,json
import numpy as np
import scipy.constants as pc
from pymatgen.io.vasp.inputs import Incar


class NBandsEsti(object):
    def __init__ (self, 
                  test_list):
        if type(test_list) is list:
            ele_t = []
            vol = []
            d_nbd = []
            nbd = []
            for ii in test_list:
                res = NBandsEsti._get_res(ii)
                ele_t.append(res['ele_temp'])
                vol.append(res['vol'])
                d_nbd.append(NBandsEsti._get_default_nbands(res))
                nbd.append(res['nbands'])
            ele_t = np.array(ele_t)
            vol = np.array(vol)
            d_nbd = np.array(d_nbd)
            nbd = np.array(nbd)
            alpha = (nbd - d_nbd) / vol / ele_t**1.5 
            self.err = np.std(alpha)
            self.pref = np.average(alpha)
            # print(np.average(alpha), np.std(alpha), self.err/self.pref)
            # print((ele_t), vol, d_nbd, nbd, alpha)
        elif type(test_list) is str:
            with open(test_list) as fp:
                self.pref = float(fp.readline())
                self.err = float(fp.readline())
        else:
            raise RuntimeError('unknown input type ' + type(test_list))

    def save(self, fname):
        with open(fname, 'w') as fp:
            fp.write(str(self.pref) + '\n')
            fp.write(str(self.err) + '\n')

    def predict(self, 
                target_dir, 
                tolerance = 0.5):
        res = NBandsEsti._get_res(target_dir)
        ele_t=(res['ele_temp'])
        vol=(res['vol'])
        d_nbd=(NBandsEsti._get_default_nbands(res))
        nbd=(res['nbands'])
        esti = (self.pref + tolerance*self.err) * ele_t**1.5 * vol + d_nbd
        return int(esti)+1

    @classmethod
    def _get_res(self, res_dir):
        res = {}
        sys = dpdata.System(os.path.join(res_dir, 'POSCAR'))
        res['natoms'] = (sys['atom_numbs'])
        res['vol'] = np.linalg.det(sys['cells'][0])
        res['nvalence'] = (self._get_potcar_nvalence(os.path.join(res_dir, 'POTCAR')))
        res['ele_temp'] = self._get_incar_ele_temp(os.path.join(res_dir, 'INCAR')) * pc.electron_volt / pc.Boltzmann
        res['nbands'] = self._get_incar_nbands(os.path.join(res_dir, 'INCAR'))
        return res

    @classmethod
    def _get_default_nbands(self, res):
        ret = 0
        for ii,jj in zip(res['natoms'], res['nvalence']):
            ret += ii * jj // 2 + ii // 2 + 2
        return ret

    @classmethod
    def _get_potcar_nvalence(self, fname):
        with open(fname) as fp:
            pot_str = fp.read().split('\n')
        head_idx = []
        for idx,ii in enumerate(pot_str):
            if ('PAW_' in ii) and ('TITEL' not in ii):
                head_idx.append(idx)
        res = []
        for ii in head_idx:
            res.append(float(pot_str[ii+1]))
        return res

    @classmethod
    def _get_incar_ele_temp(self, fname):
        incar = Incar.from_file(fname)
        return incar['SIGMA']

    @classmethod
    def _get_incar_nbands(self, fname):
        incar = Incar.from_file(fname)
        return incar.get('NBANDS')
