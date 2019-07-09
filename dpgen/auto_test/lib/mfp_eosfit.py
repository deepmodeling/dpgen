#!/usr/bin/env python3

from __future__ import division
import os
import sys
import argparse

import numpy as np
from scipy.optimize import leastsq, root, fsolve, curve_fit
from scipy.optimize import minimize
from scipy.misc import derivative
from scipy.interpolate import *
import scipy.integrate as INT
import matplotlib.pyplot as plt

kb = 1.3806488e-23  # J K^-1
kb_ev = 8.6173324e-05  # eV K^-1
h = 6.62606957e-34      # J.s
h_ev = 4.135667516e-15  # eV s
hb = 1.054571726e-34  # J.s
hb_ev = 6.58211928e-16  # eV s
mu = 1.660538921e-27  # kg
me = 9.1093821545e-31  # kg
NA = 6.02214129e+23  # number per mol

eV2GPa = 1.602176565e+2
eV2mol = 9.648455461e4  # eV/K --> J/mol/K


def __version__():
    return '1.2.5'


def get_eos_list_4p():
    eos_list_4p = [
        'murnaghan', 'birch',
        'BM4', 'mBM4', 'mBM4poly', 'rBM4',
        'rPT4', 'LOG4',
        'vinet', 'Li4p', 'universal',
        'morse', 'morse_AB', 'mie', 'mie_simple',
        'SJX_v2'
    ]
    return eos_list_4p


def get_eos_list_5p():
    eos_list_5p = [
        'BM5', 'mBM5', 'rBM5', 'mBM5poly',
        'rPT5', 'LOG5',
        'TEOS', 'SJX_5p']
    return eos_list_5p


def get_eos_list_6p():
    eos_list_6p = ['morse_6p']
    return eos_list_6p


def get_eos_list_3p():
    eos_list_3p = ['morse_3p']
    return eos_list_3p


def get_eos_list():
    LIST_ALL = get_eos_list_3p() + get_eos_list_4p() + \
        get_eos_list_5p() + get_eos_list_6p()
    return LIST_ALL


# ----------------------------------------------------------------------------------------
def res_murnaghan(pars, y, x):
    return y - murnaghan(x, pars)


def murnaghan(vol, pars):
    """
    Four-parameters murnaghan EOS.
    From PRB 28,5480 (1983)
    """
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    xx = (v0 / vol)**bp
    A = e0 - b0 * v0 / (bp - 1)
    B = b0 * vol / bp
    ee = A + B * (1 + xx / (bp - 1))

    # ee = e0+b0*vol/bp*(((v0/vol)**bp)/(bp-1)+1) -b0*v0/(bp-1)

    return ee

# ----------------------------------------------------------------------------------------


def res_birch(pars, y, x):
    return y - birch(x, pars)


def birch(v, parameters):
    """
    From Intermetallic compounds: Principles and Practice, Vol. I: Princples
    Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos
    paper downloaded from Web

    case where n=0
    """
    e0 = parameters[0]
    b0 = parameters[1]
    bp = parameters[2]
    v0 = parameters[3]

    e = (e0 + 9.0 / 8.0 * b0 * v0 * ((v0 / v) ** (2.0 / 3.0) - 1.0) ** 2
         + 9.0 / 16.0 * b0 * v0 * (bp - 4.) * ((v0 / v) ** (2.0 / 3.0) - 1.0) ** 3)

    return e


# ----------------------------------------------------------------------------------------
def res_mBM4(pars, y, x):
    return y - mBM4(x, pars)


def calc_props_mBM4(pars):
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    bpp = (-74 + 45 * bp - 9 * bp * bp) / (9 * b0)
    props = [e0, b0, bp, v0, bpp]
    return props


def mBM4(vol, pars):
    """
    Birch-Murnaghan 4 pars equation from PRB 70, 224107, 3-order BM
    """
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    a = e0 + 9 * b0 * v0 * (4 - bp) / 2.
    b = -9 * b0 * v0 ** (4. / 3) * (11 - 3 * bp) / 2.
    c = 9 * b0 * v0 ** (5. / 3) * (10 - 3 * bp) / 2.
    d = -9 * b0 * v0 ** 2 * (3 - bp) / 2.

    n = 1  # 1 as mBM, 2 as BM
    VV = np.power(vol, -n / 3)
    ee = a + b * VV + c * VV ** 2 + d * VV ** 3

    return ee


def res_mBM5(pars, y, x):
    return y - mBM5(x, pars)


def mBM5(vol, pars):
    """
    modified BM5 EOS, Shang SL comput mater sci, 2010: 1040-1048
    """
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]
    b2p = pars[4]

    '''
    # copy from ShunLi's matlab scripts.
    a = (8 * e0 + 3 * b0 * (122 + 9 * b0 * b2p - 57 * bp + 9 * bp * bp) * v0) / 8
    b = (-3 * b0 * (107 + 9 * b0 * b2p - 54 * bp + 9 * bp * bp) * v0 ** (4 / 3)) / 2
    c = (9 * b0 * (94 + 9 * b0 * b2p - 51 * bp + 9 * bp * bp) * v0 ** (5 / 3)) / 4
    d = (-3 * b0 * (83 + 9 * b0 * b2p - 48 * bp + 9 * bp * bp) * v0 ** 2) / 2
    e = (3 * b0 * (74 + 9 * b0 * b2p - 45 * bp + 9 * bp * bp) * v0 ** (7 / 3)) / 8
    '''
    # retype according to formula in the article
    a = e0 + 3 * b0 * v0 * (122 + 9 * b0 * b2p - 57 * bp + 9 * bp * bp) / 8
    b = -3 * b0 * v0**(4 / 3) * (107 + 9 * b0 *
                                 b2p - 54 * bp + 9 * bp * bp) / 2
    c = 9 * b0 * v0**(5 / 3) * (94 + 9 * b0 * b2p - 51 * bp + 9 * bp * bp) / 4
    d = -3 * b0 * v0**2 * (83 + 9 * b0 * b2p - 48 * bp + 9 * bp * bp) / 2
    e = 3 * b0 * v0**(7 / 3) * (74 + 9 * b0 * b2p - 45 * bp + 9 * bp * bp) / 8

    VV = np.power(vol, -1 / 3)
    ee = a + b * VV + c * VV ** 2 + d * VV ** 3 + e * VV ** 4

    return ee


# ----------------------------------------------------------------------------------------
def res_mBM4poly(pars, y, x):
    return y - mBM4poly(x, pars)


def mBM4poly(vol, parameters):
    """
    modified BM5 EOS, Shang SL comput mater sci, 2010: 1040-1048, original expressions.
    """
    a = parameters[0]
    b = parameters[1]
    c = parameters[2]
    d = parameters[3]
    e = 0

    n = 1  # 1 as mBM, 2 as BM
    VV = np.power(vol, -1 / 3)
    E = a + b * VV + c * VV ** 2 + d * VV ** 3 + e * VV ** 4
    return E


def calc_v0_mBM4poly(x, pars):
    a = pars[0]
    b = pars[1]
    c = pars[2]
    d = pars[3]
    e = 0

    f = ((4 * e) / (3 * x ** (7 / 3)) + d / x ** 2 + (2 * c) /
         (3 * x ** (5 / 3)) + b / (3 * x ** (4 / 3))) * eV2GPa
    return f


def calc_props_mBM4poly(pars):
    a = pars[0]
    b = pars[1]
    c = pars[2]
    d = pars[3]
    e = 0

    v0 = 4 * c ** 3 - 9 * b * c * d + \
        np.sqrt((c ** 2 - 3 * b * d) * (4 * c ** 2 - 3 * b * d) ** 2)
    v0 = -v0 / b ** 3
    b0 = ((28 * e) / (9 * v0 ** (10 / 3)) +
          (2 * d) / v0 ** 3 + (10 * c) / (9 * v0 ** (8 / 3)) +
          (4 * b) / (9 * v0 ** (7 / 3))) * v0
    bp = (98 * e + 54 * d * v0 ** (1 / 3) + 25 * c * v0 ** (2 / 3) + 8 * b * v0) / \
        (42 * e + 27 * d * v0 ** (1 / 3) + 15 * c * v0 ** (2 / 3) + 6 * b * v0)
    b2p = (v0**(8 / 3) * (9 * d * (14 * e + 5 * c * v0**(2 / 3) + 8 * b * v0) +
                          2 * v0**(1 / 3) * (126 * b * e * v0 ** (1 / 3) + 5 * c * (28 * e + b * v0)))) / \
        (2 * (14 * e + 9 * d * v0**(1 / 3) + 5 * c * v0**(2 / 3) + 2 * b * v0)**3)
    e0 = mBM4poly(v0, pars)

    props = [e0, b0, bp, v0, b2p]
    return props

#---------------------------------


def res_mBM5poly(pars, y, x):
    return y - mBM5poly(x, pars)


def mBM5poly(vol, pars):
    """
    modified BM5 EOS, Shang SL comput mater sci, 2010: 1040-1048, original expressions.
    """
    a = pars[0]
    b = pars[1]
    c = pars[2]
    d = pars[3]
    e = pars[4]

    n = 1  # 1 as mBM, 2 as BM
    VV = np.power(vol, -1 / 3)
    E = a + b * VV + c * VV ** 2 + d * VV ** 3 + e * VV ** 4
    return E


def calc_v0_mBM5poly(x, pars):
    a = pars[0]
    b = pars[1]
    c = pars[2]
    d = pars[3]
    e = pars[4]

    f = ((4 * e) / (3 * x ** (7 / 3)) + d / x ** 2 + (2 * c) /
         (3 * x ** (5 / 3)) + b / (3 * x ** (4 / 3))) * eV2GPa
    return f


def calc_props_mBM5poly(pars):
    a = pars[0]
    b = pars[1]
    c = pars[2]
    d = pars[3]
    e = pars[4]
    n = 1

    # guess v0
    vtest = np.linspace(1, 100, 101)
    etest = mBM5poly(vtest, pars)
    vg = vtest[etest.argmin()]

    # sol = root(calc_v0_mBM5poly, vg, args=pars)
    # v0 = sol.x[0]
    sol = fsolve(calc_v0_mBM5poly, vg, args=pars)
    v0 = sol

    b0 = ((28 * e) / (9 * v0 ** (10 / 3)) +
          (2 * d) / v0 ** 3 + (10 * c) / (9 * v0 ** (8 / 3)) +
          (4 * b) / (9 * v0 ** (7 / 3))) * v0
    bp = (98 * e + 54 * d * v0 ** (1 / 3) + 25 * c * v0 ** (2 / 3) + 8 * b * v0) / \
        (42 * e + 27 * d * v0 ** (1 / 3) + 15 * c * v0 ** (2 / 3) + 6 * b * v0)
    b2p = (v0**(8 / 3) * (9 * d * (14 * e + 5 * c * v0**(2 / 3) + 8 * b * v0) +
                          2 * v0**(1 / 3) * (126 * b * e * v0 ** (1 / 3) + 5 * c * (28 * e + b * v0)))) / \
        (2 * (14 * e + 9 * d * v0**(1 / 3) + 5 * c * v0**(2 / 3) + 2 * b * v0)**3)

    '''
    dEdV = -b / 3 * v0**(-4 / 3) - 2 / 3 * c * \
        v0**(-5 / 3) - d * v0**(-2) - 4 / 3 * e * v0**(-7 / 3)
    dEdV2 = 4 / 9 * b * v0**(-7 / 3) + 10 / 9 * c * \
        v0**(-8 / 3) + 2 * d * v0**(-3) + 28 / 9 * e * v0**(-10 / 3)
    dEdV3 = -28 / 27 * b * \
        v0**(-10 / 3) - 80 / 27 * c * v0**(-11 / 3) - 6 * \
        d * v0**(-4) - 280 / 27 * e * v0**(-13 / 3)
    dEdV4 = 280 / 81 * b * v0**(-13 / 3) + 880 / 81 * c * \
        v0**(-14 / 3) + 24 * d * v0**(-5) + 3640 / 81 * e * v0**(-16 / 3)

    P = -v0 * dEdV
    dBdV = dEdV2 + v0 * dEdV3
    dPdV = -dEdV - v0 * dEdV2
    dBdV2 = dEdV3 + dEdV3 + v0 * dEdV4
    dPdV2 = -dEdV2 - dEdV2 - v0 * dEdV3

    b0 = v0 * dEdV2
    bp = dBdV / dPdV
    b2p = (dBdV2 * dPdV - dPdV2 * dBdV) / dPdV**3
    '''
    e0 = mBM5poly(v0, pars)
    props = [e0, b0, bp, v0, b2p]
    return props


# ----------------------------------------------------------------------------------------


def res_BM4(pars, y, x):
    return y - BM4(x, pars)


def calc_props_BM4(pars):
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]
    bpp = (-143 + 63 * bp - 9 * bp * bp) / (9 * b0)
    props = [e0, b0, bp, v0, bpp]
    return props


def BM4(vol, pars):
    """
    Birch-Murnaghan 4 pars equation from PRB 70, 224107, 3-order
    """
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    eta = (v0 / vol) ** (1. / 3.)
    e = e0 + 9. * b0 * v0 / 16 * \
        (eta ** 2 - 1) ** 2 * (6 + bp * (eta ** 2 - 1.) - 4. * eta ** 2)
    return e


# ----------------------------------------------------------------------------------------
def res_BM5(pars, y, x):
    return y - BM5(x, pars)


def BM5(vol, pars):
    """
    Birch-Murnaghan 5 pars equation from PRB 70, 224107, 4-Order
    """
    e0 = pars[0]
    b0 = pars[1]
    b0p = pars[2]
    v0 = pars[3]
    b0pp = pars[4]

    t1 = (v0 / vol) ** (1. / 3.)
    t2 = t1 ** 2
    t3 = t2 - 1.
    t4 = t3 ** 2 / 4.
    t5 = b0p ** 2

    ee = e0 + 3. / 8. * b0 * v0 * t4 * (9. * t4 * b0 * b0pp +
                                        9. * t4 * t5 - 63. * t4 * b0p + 143. * t4 + 6. * b0p * t3 - 24. * t2 + 36.)

    return ee


def rBM4(vol, pars):
    '''
    Implementions as Alberto Otero-de-la-Roza, i.e. rBM4 is used here
    Comput Physics Comm, 2011, 182: 1708-1720
    '''
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    x = v0 / vol
    f = 0.5 * (x**(2. / 3) - 1)
    E = e0 + 4.5 * v0 * b0 * f**2 * (1 + (bp - 4) * f)
    return E


def res_rBM4(pars, y, x):
    res = y - rBM4(x, pars)
    return res


def rBM4_pv(vol, pars):
    '''
    Implementions as Alberto Otero-de-la-Roza, i.e. rBM4 is used here
    Comput Physics Comm, 2011, 182: 1708-1720
    Fit for V-P relations
    '''
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    x = v0 / vol
    f = 0.5 * (x**(2. / 3) - 1)
    P = 1.5 * b0 * (2 * f + 1)**(2.5) * (2 + 3 * (bp - 4) * f)
    return P


def res_rBM4_pv(par, y, x):
    res = y - rBM4_pv(x, pars)
    return res


def rBM5(vol, pars):
    '''
    Implementions as Alberto Otero-de-la-Roza, i.e. rBM5 is used here
    Comput Physics Comm, 2011, 182: 1708-1720
    '''
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]
    bpp = pars[4]

    x = v0 / vol
    f = 0.5 * (x**(2. / 3) - 1)
    H = b0 * bpp + b0 * b0

    E = e0 + 3 / 8 * v0 * b0 * f**2 * \
        ((9 * H - 63 * bp + 143) * f**2 + 12 * (bp - 4) * f + 12)
    return E


def res_rBM5(pars, y, x):
    res = y - rBM5(x, pars)
    return res


def rBM5_pv(vol, pars):
    '''
    Implementions as Alberto Otero-de-la-Roza, i.e. rBM5 is used here
    Comput Physics Comm, 2011, 182: 1708-1720
    Fit for V-P relations
    '''
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]
    bpp = pars[4]

    x = v0 / vol
    f = 0.5 * (x**(2. / 3) - 1)
    H = b0 * bpp + b0 * b0

    P = 0.5 * b0 * (2 * f + 1)**(2.5) * \
        ((9 * H - 63 * bp + 143) * f**2 + 9 * (bp - 4) * f + 6)
    return P


def res_rBM5_pv(par, y, x):
    res = y - rBM5_pv(x, pars)
    return res


# ----------------------------------------------------------------------------------------


def res_universal(pars, y, x):
    return y - universal(x, pars)


def universal(vol, parameters):
    """
    Universal equation of state(Vinet P et al., J. Phys.: Condens. Matter 1, p1941 (1989))
    """
    e0 = parameters[0]
    b0 = parameters[1]
    bp = parameters[2]
    v0 = parameters[3]

    t1 = b0 * v0
    t2 = bp - 1.
    t3 = (vol / v0) ** (1. / 3.)
    t4 = np.exp(-3. / 2. * t2 * (-1. + t3))
    t5 = t2 ** 2
    t6 = 1. / t5
    e = e0 - 2. * t1 * t4 * \
        (3. * t3 * bp - 3. * t3 + 5. - 3. * bp) * t6 + 4. * t1 * t6

    return e


# ----------------------------------------------------------------------------------------
def res_LOG4(pars, y, x):
    return y - LOG4(x, pars)


def LOG4(vol, pars):
    '''
    Natrual strain (Poirier-Tarantola)EOS with 4 paramters
    Seems only work in near-equillibrium range.
    '''
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    t1 = b0 * v0
    t2 = np.log(v0 / vol)
    t3 = t2 ** 2
    t4 = t3 * t2
    ee = e0 + t1 * t3 / 2. + t1 * t4 * bp / 6. - t1 * t4 / 3.
    '''
    # write follows ShunLi's
    xx = np.log(v0)
    a = e0 + b0 * v0 * (3 * xx**2 + (bp - 2) * xx**3) / 6
    b = -b0 * v0 * (2 * xx + (bp - 2) * xx**2) / 2
    c = b0 * v0 * (1 + (bp - 2) * xx) / 2
    d = -b0 * v0 * (bp - 2) / 6
    VV = np.log(vol)
    ee = a + b * VV + c * VV**2 + d * VV**3
    '''
    return ee


def calc_props_LOG4(pars):
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    bpp = (-3 + 3 * bp - bp * bp) / b0
    props = [e0, b0, bp, v0, bpp]
    return props

# ----------------------------------------------------------------------------------------


def rPT4(vol, pars):
    '''
    Natrual strain EOS with 4 paramters
    Seems only work in near-equillibrium range.
    Implementions as Alberto Otero-de-la-Roza, i.e. rPT4 is used here
    Comput Physics Comm, 2011, 182: 1708-1720,
    in their article, labeled as PT3 (3-order), however, we mention it as
    rPT4 for 4-parameters EOS.
    '''
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    x = vol / v0
    fn = 1 / 3 * np.log(x)
    E = e0 + 4.5 * b0 * v0 * fn**2 * (-(bp - 2) * fn + 1)
    return E


def res_rPT4(pars, y, x):
    res = y - rPT4(x, pars)
    return res


def rPT4_pv(vol, pars):
    '''
    Natrual strain (Poirier-Tarantola)EOS with 4 paramters
    Seems only work in near-equillibrium range.
    Implementions as Alberto Otero-de-la-Roza, i.e. rPT4 is used here
    Comput Physics Comm, 2011, 182: 1708-1720,
    in their article, labeled as PT3 (3-order), however, we mention it as
    rPT4 for 4-parameters EOS.
    '''
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    x = (vol / v0)
    fn = 1 / 3 * np.log(x)
    P = -1.5 * b0 * fn * np.exp(-3 * fn) * (-3 * (bp - 2) * fn + 1)
    return P


def res_rPT4_pv(pars, y, x):
    res = y - rPT4_pv(x, pars)
    return res


# ----------------------------------------------------------------------------------------
def res_LOG5(pars, y, x):
    return y - LOG5(x, pars)


def LOG5(vol, parameters):
    '''
    Natrual strain (Poirier-Tarantola)EOS with 5 paramters
    '''
    e0 = parameters[0]
    b0 = parameters[1]
    b0p = parameters[2]
    v0 = parameters[3]
    b0pp = parameters[4]

    t1 = b0 * v0
    t2 = np.log(v0 / vol)
    t3 = t2 ** 2
    t4 = t3 ** 2
    t5 = b0 ** 2
    t6 = b0p ** 2
    t7 = t3 * t2
    e = e0 + t1 * t4 / 8. + t5 * v0 * t4 * b0pp / 24. - t1 * t4 * b0p / 8. + \
        t1 * t4 * t6 / 24. + t1 * t7 * b0p / 6. - t1 * t7 / 3. + t1 * t3 / 2.

    return e


def rPT5(vol, pars):
    '''
    Natrual strain EOS with 4 paramters
    Seems only work in near-equillibrium range.
    Implementions as Alberto Otero-de-la-Roza, i.e. rPT5 is used here
    Comput Physics Comm, 2011, 182: 1708-1720,
    in their article, labeled as PT3 (3-order), however, we mention it as
    rPT5 for 4-parameters EOS.
    '''
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]
    bpp = pars[4]

    x = (vol / v0)
    fn = 1 / 3 * np.log(x)
    H = b0 * bpp + bp * bp

    E = e0 + 9 / 8 * b0 * v0 * fn**2 * \
        (-3 * (H + 3 * bp - 3) * fn**2 - 4 * (bp - 2) * fn + 4)
    return E


def res_rPT5(pars, y, x):
    res = y - rPT5(x, pars)
    return res


def rPT5_pv(vol, pars):
    '''
    Natrual strain (Poirier-Tarantola)EOS with 5 paramters
    Implementions as Alberto Otero-de-la-Roza, i.e. rPT5 is used here
    Comput Physics Comm, 2011, 182: 1708-1720,
    in their article, labeled as PT3 (3-order), however, we mention it as
    rPT5 for 4-parameters EOS.
    '''
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]
    bpp = pars[4]

    x = (vol / v0)
    fn = 1 / 3 * np.log(x)
    H = b0 * bpp + bp * bp

    P = -1.5 * b0 * fn * \
        np.exp(-3 * fn) * (-3 * (H + 3 * bp - 3)
                           * fn**2 - 3 * (bp - 2) * fn + 2)
    return P


def res_rPT5_pv(pars, y, x):
    res = y - rPT5_pv(x, pars)
    return res


# ----------------------------------------------------------------------------------------
def res_vinet(pars, y, x):
    res = y - vinet(x, pars)
    return res


def calc_props_vinet(pars):
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    bpp = (19 - 18 * bp - 9 * bp * bp) / (36 * b0)
    props = [e0, b0, bp, v0, bpp]
    return props


def vinet(vol, pars):
    """
    Vinet equation from PRB 70, 224107
    Following, Shang Shunli et al., comput mater sci, 2010: 1040-1048, original expressions.
    """
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    eta = (vol / v0) ** (1 / 3)
    Y = 1.5 * (bp - 1) * (1 - eta)
    Z = 4 * b0 * v0 / (bp - 1)**2

    ee = e0 + Z - Z * (1 - Y) * np.exp(Y)
    return ee


def vinet_pv(vol, pars):
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    eta = (vol / v0) ** (1 / 3)
    P = 3 * b0 * (1 - eta / eta**2) * np.exp(-1.5 * (bp - 1) * (eta - 1))
    return P


def res_vinet_pv(par, y, x):
    res = y - vinet(x, pars)
    return res


# ----------------------------------------------------------------------------------------
def Li4p(V, parameters):
    ''' Li JH, APL, 87, 194111 (2005) '''
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]

    # In fact, it is just a small modified version of Rose or Vinet EOS.
    x = (V / V0) ** (1. / 3.)
    eta = np.sqrt(-9 * B0 * V0 / E0)
    astar = eta * (x - 1.)
    delta = (BP - 1.) / (2 * eta) - 1. / 3.

    E = E0 * (1 + astar + delta * astar ** 3) * np.exp(-astar)

    return E


def res_Li4p(p, y, x):
    return y - Li4p(x, p)


# ----------------------------------------------------------------------------------------
def morse(v, pars):
    ''' Reproduce from ShunliShang's matlab script. '''
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]

    # usually used in high pressure range, but maybe failed in large
    # expansion( > 2*v0) range.
    a = e0 + (9 * b0 * v0) / (2 * (-1 + bp) * (-1 + bp))
    b = (-9 * b0 * np.exp(-1 + bp) * v0) / (-1 + bp) / (-1 + bp)
    c = (9 * b0 * np.exp(-2 + 2 * bp) * v0) / (2 * (-1 + bp) * (-1 + bp))
    d = (1 - bp) / (v0 ** (1. / 3))

    ee = a + b * np.exp(d * np.power(v, 1. / 3)) + \
        c * np.exp(2 * d * np.power(v, 1. / 3))
    return ee


def calc_props_morse(pars):
    e0 = pars[0]
    b0 = pars[1]
    bp = pars[2]
    v0 = pars[3]
    bpp = (5 - 5 * bp - 2 * bp * bp) / (9 * b0)

    props = [e0, b0, bp, v0, bpp]
    return props


def res_morse(p, en, volume):
    res = en - morse(volume, p)
    return res


def morse_AB(volume, p):
    '''
    morse_AB EOS formula from Song's FVT souces
    '''
    # p0 = [e0, b0, bp, v0, bpp]
    E0 = p[0]
    A = p[1]
    B = p[2]
    V0 = p[3]
    xx = (volume / V0) ** (1 / 3)
    E = 1.5 * V0 * A * (np.exp(B * (1. - 2.0 * xx)) -
                        2.0 * np.exp(-B * xx) + np.exp(-B)) + E0
    return E


def res_morse_AB(p, en, volume):
    res = en - morse_AB(volume, p)
    return res


def morse_3p(volume, p):
    '''
    morse_AB EOS formula from Song's FVT souces
    A= 0.5*B
    '''
    # p0 = [e0, b0, bp, v0, bpp]
    E0 = p[0]
    A = p[1]
    V0 = p[2]
    B = 0.5 * A
    xx = (volume / V0) ** (1 / 3)
    E = 1.5 * V0 * A * (np.exp(B * (1. - 2.0 * xx)) -
                        2.0 * np.exp(-B * xx) + np.exp(-B)) + E0
    return E


def res_morse_3p(p, en, volume):
    res = en - morse_3p(volume, p)
    return res


def morse_6p(vol, par):
    '''
    Generalized Morse EOS proposed by Qin, see:
    Qin et al. Phys Rev B, 2008, 78, 214108.
    Qin et al. Phys Rev B, 2008, 77, 220103(R).
    '''
    # p0 = [e0, b0, bp, v0, bpp]
    e0 = par[0]
    b0 = par[1]
    p = abs(par[2])
    v0 = par[3]
    q = abs(par[4])
    m = abs(par[5])
    n = abs(par[6])

    C = (p + m) * (q - n) * ((p + m) - (q - n)) + (p * n + q * m)
    A = 9 * b0 * v0 * (q - n) / C
    B = 9 * b0 * v0 * (p + m) / C

    x = (vol / v0) ** (1 / 3)
    ee = A * np.exp(-p * (x - 1)) / (x ** m) - \
        B * (x ** n) * np.exp(-q * (x - 1)) + e0 - (A - B)
    return ee


def calc_props_morse_6p(par):
    e0 = par[0]
    b0 = par[1]
    p = abs(par[2])
    v0 = par[3]
    q = abs(par[4])
    m = abs(par[5])
    n = abs(par[6])

    bp = (p + q) / 3 + 1
    bpp = ((7 - p * q) / 9 - bp) / b0
    props = [e0, b0, bp, v0, bpp, m, n]
    return props


def res_morse_6p(p, en, volume):
    res = en - morse_6p(volume, p)
    return res


# ----------------------------------------------------------------------------------------
def mie(v, p):
    '''
    Mie model for song's FVT
    '''
    # p0 = [e0, b0, bp, v0, bpp]
    E0 = p[0]
    m = p[1]
    n = p[2]
    V0 = p[3]

    xx = (V0 / v) ** (1 / 3)
    E = E0 / (n - m) * (n * (xx ** m) - m * (xx ** n))

    return E


def res_mie(p, e, v):
    res = e - mie(v, p)
    return res


def mie_simple(v, p):
    '''
    Mie_simple model for song's FVT
    '''
    # p0 = [e0, b0, bp, v0, bpp]
    E0 = p[0]
    m = 4
    n = p[2]
    V0 = p[3]

    xx = (V0 / v) ** (1 / 3)
    E = E0 / (n - m) * (n * (xx ** m) - m * (xx ** n))

    return E


def res_mie_simple(p, e, v):
    res = e - mie_simple(v, p)
    return res


# ----------------------------------------------------------------------------------------
def TEOS(v, par):
    '''
    Holland, et al, Journal of Metamorphic Geology, 2011, 29(3): 333-383
    Modified Tait equation of Huang & Chow
    '''
    e0 = par[0]
    b0 = par[1]
    bp = par[2]
    v0 = par[3]
    bpp = par[4]

    a = (1 + bp) / (1 + bp + b0 * bpp)
    b = bp / b0 - bpp / (1 + bp)
    c = (1 + bp + b0 * bpp) / (bp ** 2 + bp - b0 * bpp)

    t1 = (v / v0) + a - 1
    P = 1 / b * ((t1 / a) ** (-1 / c) - 1)
    # ?? the following energy relation is nonsenses.
    E = -P * np.log(v0) + e0
    return E


def res_TEOS(p, e, v):
    res = e - TEOS(v, p)
    return res


# ----------------------------------------------------------------------------------------
def SJX_v2(vol, par):
    '''
    Sun Jiuxun, et al. J phys Chem Solids, 2005, 66: 773-782.
    They said it is satified for the limiting condition at high pressure.
    '''
    e0 = par[0]
    b0 = par[1]
    bp = par[2]
    v0 = par[3]

    Y = (v0 / vol)**(1 / 3) - 1
    alpha = 1 / 4 * (3 * bp - 1)
    beta = 1 - 1 / alpha
    a = alpha
    b = beta
    ee = e0 + 3 / 10 * a**4 * b0 * v0 * (2 * (b**5 + 1 / a**5) * (1 - (Y + 1)**(-3))
                                         - 15 * b**4 * (1 - (Y + 1)**(-2))
                                         - 60 * b**3 * (1 - (Y + 1)**(-1))
                                         + 60 * b**2 * np.log(Y + 1) + 3 * Y * (Y + 2) - 30 * b * Y)
    return ee


def res_SJX_v2(p, e, v):
    res = e - SJX_v2(v, p)
    return res


def SJX_5p(vol, par):
    '''
    SJX_5p's five parameters EOS, Physica B: Condens Mater, 2011, 406: 1276-1282
    '''
    e0 = par[0]
    a = par[1]
    b = par[2]
    v0 = par[3]
    n = par[4]

    X = (vol / v0) ** (1 / 3)
    C1 = n * b * np.exp(a * (1 - X) + b * (1 - X ** n))
    C2 = (-n * b - a) * np.exp(b * (1 - X ** n))
    ee = e0 / a * (C1 + C2)
    return ee


def res_SJX_5p(p, e, v):
    res = e - SJX_5p(v, p)
    return res


def calc_props_SJX_5p(par):
    e0 = par[0]
    a = par[1]
    b = par[2]
    v0 = par[3]
    n = par[4]

    b0 = n * b * e0 / (9 * v0) * (a + n * b + n - 1)
    Q = (a + n * b) * (a + 2 * n * b) - (n - 1) * (n - 2)
    bp = 1 + (Q / 3) / (a + n * b + n - 1)

    props = [e0, b0, bp, v0, n]
    return props

# ----------------------------------------------------------------------------------------
# some utility functions


def read_ve(fin):
    if not os.path.exists(fin):
        print('Could not find input file: [%s]' % fin)
        os.sys.exit(-1)
    lines = open(fin).readlines()
    nline = len(lines)
    vol = []
    en = []
    for i in range(nline):
        tmp = lines[i].split()
        if len(tmp) >= 2:
            # now, we can read ve.dat & vec.dat
            # read the first two numbers.
            v = float(tmp[0])
            e = float(tmp[1])
            vol.append(v)
            en.append(e)
        else:
            pass
    return [vol, en]


def read_vlp(fin, fstart, fend):
    if not os.path.exists(fin):
        print('>> Could not find input file: [%s]' % fin)
        os.sys.exit(-1)
    lines = open(fin).readlines()
    nline = len(lines)
    print("\n** Totally, there are %d data points reading..." % nline)
    vol = []
    cella = []
    cellb = []
    cellc = []
    cellba = []
    cellca = []
    for i in range(nline):
        tmp = lines[i].split()
        if len(tmp) >= 2:
            # read the DATA information
            v = float(tmp[0])
            if len(tmp) == 6:
                # for vlp.dat
                a = float(tmp[1])
                b = float(tmp[2])
                c = float(tmp[3])
                ba = float(tmp[4])
                ca = float(tmp[5])
            elif len(tmp) == 7:
                # for velp.dat
                e = float(tmp[1])
                a = float(tmp[2])
                b = float(tmp[3])
                c = float(tmp[4])
                ba = float(tmp[5])
                ca = float(tmp[6])
            vol.append(v)
            cella.append(a)
            cellb.append(b)
            cellc.append(c)
            cellba.append(ba)
            cellca.append(ca)
    print("\n** Vmin = %f, Vmax = %f" % (min(vol), max(vol)))

    # some special conditions
    if (fstart <= 0):
        print('\n** Data range input parameters must be positive values!')
        os.sys.exit(-1)
    if fstart > fend:
        if fend != -1:
            tmp = fstart
            fstart = fend
            fend = tmp
    if (fend > nline):
        print(
            '\n** EoSfit fit range exceed available data numbers, Reset it to be %d now.' % nline)
        fend = nline
    if (fstart > nline):
        print('EoSfit fit range exceed available data numbers, Reset it to be 1: %d now.' % nline)
        fstart = 1
        fend = -1

    fstart = fstart - 1

    # get data in the predefined range.
    fstart = fstart - 1
    if fend != -1:
        vol = vol[fstart:fend]
        cella = cella[fstart:fend]
        cellb = cellb[fstart:fend]
        cellc = cellc[fstart:fend]
        cellba = cellba[fstart:fend]
        cellca = cellca[fstart:fend]
    else:
        vol = vol[fstart:]
        cella = cella[fstart:]
        cellb = cellb[fstart:]
        cellc = cellc[fstart:]
        cellba = cellba[fstart:]
        cellca = cellca[fstart:]
    return [vol, cella, cellb, cellc, cellba, cellca]


def read_velp(fin, fstart, fend):
    if not os.path.exists(fin):
        print('>> Could not find input file: [%s]' % fin)
        os.sys.exit(-1)
    lines = open(fin).readlines()
    nline = len(lines)
    print("\n** Totally, there are %d data points reading..." % nline)
    vol = []
    eng = []
    cella = []
    cellb = []
    cellc = []
    cellba = []
    cellca = []
    for i in range(nline):
        tmp = lines[i].split()
        if len(tmp) >= 2:
            # read the DATA information
            v = float(tmp[0])
            if len(tmp) == 7:
                # for velp.dat
                e = float(tmp[1])
                a = float(tmp[2])
                b = float(tmp[3])
                c = float(tmp[4])
                ba = float(tmp[5])
                ca = float(tmp[6])
            vol.append(v)
            eng.append(e)
            cella.append(a)
            cellb.append(b)
            cellc.append(c)
            cellba.append(ba)
            cellca.append(ca)
    print("\n** Vmin = %f, Vmax = %f" % (min(vol), max(vol)))

    # some special conditions
    if (fstart <= 0):
        print('\n** Data range input parameters must be positive values!')
        os.sys.exit(-1)
    if fstart > fend:
        if fend != -1:
            tmp = fstart
            fstart = fend
            fend = tmp
    if (fend > nline):
        print(
            '\n** EoSfit fit range exceed available data numbers, Reset it to be %d now.' % nline)
        fend = -1
    if (fstart > nline):
        print('EoSfit fit range exceed available data numbers, Reset it to be 1: %d now.' % nline)
        fstart = 1
        fend = -1

    # get data in the predefined range.
    fstart = fstart - 1
    if fend != -1:
        vol = vol[fstart:fend]
        eng = eng[fstart:fend]
        cella = cella[fstart:fend]
        cellb = cellb[fstart:fend]
        cellc = cellc[fstart:fend]
        cellba = cellba[fstart:fend]
        cellca = cellca[fstart:fend]
    else:
        vol = vol[fstart:]
        eng = eng[fstart:]
        cella = cella[fstart:]
        cellb = cellb[fstart:]
        cellc = cellc[fstart:]
        cellba = cellba[fstart:]
        cellca = cellca[fstart:]

    return [vol, eng, cella, cellb, cellc, cellba, cellca]


def init_guess(fin):
    v, e = read_ve(fin)
    a, b, c = np.polyfit(v, e, 2)  # this comes from pylab
    # initial guesses.
    v0 = np.abs(-b / (2 * a))
    e0 = a * v0 ** 2 + b * v0 + c
    b0 = 2 * a * v0
    bp = 3.0
    bpp = 1 * eV2GPa
    x0 = [e0, b0, bp, v0, bpp]
    return x0


def repro_ve(func, vol_i, p):
    ndata = len(vol_i)
    eni = np.zeros(ndata)
    for i in range(ndata):
        eni[i] = eval(func)(vol_i[i], p)
    return eni


def repro_vp(func, vol_i, pars):
    dv = 1e-5
    efunc = eval(func)
    ndata = len(vol_i)
    Press = np.zeros(ndata)
    for i in range(ndata):
        v = vol_i[i]
        P = (efunc(v + 0.5 * dv, pars) - efunc(v - 0.5 * dv, pars)) / dv
        Press[i] = -P * eV2GPa
    return Press


def ext_vec(func, fin, p0, fs, fe, vols=None, vole=None, ndata=101, refit=0, show_fig=False):
    '''
    extrapolate the data points for E-V based on the fitted parameters in small or
    very large volume range.
    '''
    # read fitted-parameters
    # pars = np.loadtxt(fin, dtype=float)
    fout = 'EoSfit.out'
    pars = lsqfit_eos(func, fin, p0, fs, fe, fout, refit=refit)

    vol, eng, cella, cellb, cellc, cellba, cellca = read_velp(fin, fs, fe)
    sca = ext_splint(vol, cellca)

    # define extrapolate range
    print("\n** Vext_start = %f, Vext_end = %f, N_ext = %d" %
          (vols, vole, ndata))
    vol_ext = np.linspace(vols, vole, ndata)
    en_ext = eval(func)(vol_ext, pars)

    # en_ext = np.zeros(ndata)
    fout = 'ext_ve_' + func + '.dat'
    fw = open(fout, 'w+')
    fw.write("%d\n" % ndata)
    for i in range(ndata):
        vx = vol_ext[i]
        ex = en_ext[i]
        if vx <= min(vol):
            cax = cellca[0]
        elif vx >= max(vol):
            cax = cellca[-1]
        else:
            cax = sca(vx)
        fw.write("%f\t%f\t%f\n" % (vx, ex, cax))
        fw.flush()
    fw.close()

    # plot the results
    vol, en = read_ve(fin)
    plt.plot(vol, en, 'o-', vol_ext, en_ext, 'rd')
    plt.legend(['dft_calc', func + '_ext'], loc='best')
    plt.savefig('ext_ve_' + func + '.png')
    if show_fig:
    	plt.show()
    plt.close()
    print("\n>> Storing the extrapolate results in %s\n" % fout)
    print("\n>> DONE!")

    return


def ext_splint(xp, yp, order=3, method='unispl'):
    if method == 'interp1d':
        # old-implement of 1-D interpolation
        # could not used in extrapolate
        SPLINT = interp1d
        return SPLINT(xp, yp, order, bounds_error=False)
    elif method == 'piecepoly':
        SPLINT = PiecewisePolynomial
        return SPLINT(xp, yp, order)
    else:
        if method == 'unispl':
            # 1-D smoothing spline fit to a given set of data
            SPLINT = UnivariateSpline
        else:
            SPLINT = LSQUnivariateSpline
        return SPLINT(xp, yp, k=order)


def ext_velp(fin, fstart, fend, vols, vole, ndata, order=3, method='unispl', fout='ext_velp.dat', show_fig=False):
    '''
    extrapolate the lattice parameters based on input data
    '''
    # read file
    vol, eng, cella, cellb, cellc, cellba, cellca = read_velp(
        fin, fstart, fend)

    # define extrapolate range
    print("\n** Vext_start = %f, Vext_end = %f, N_ext = %d" %
          (vols, vole, ndata))
    vv = np.linspace(vols, vole, ndata)

    # spline order = 3 by default
    se = ext_splint(vol, eng, order, method)
    ee = se(vv)
    sa = ext_splint(vol, cella, order, method)
    cellaa = sa(vv)
    sb = ext_splint(vol, cellb, order, method)
    cellbb = sb(vv)
    sc = ext_splint(vol, cellc, order, method)
    cellcc = sc(vv)
    sba = ext_splint(vol, cellba, order, method)
    cellbaba = sba(vv)
    sca = ext_splint(vol, cellca, order, method)
    cellcaca = sca(vv)
    cellca_cal = cellcc / cellaa

    # plot the extrapolate results
    nfigure = 6
    lp_ext = [ee, cellaa, cellbb, cellcc, cellbaba, cellcaca, cellca_cal]
    lp_ori = [eng, cella, cellb, cellc, cellba,
              cellca, np.array(cellc) / np.array(cella)]
    lp_ylabel = ['E(eV)', 'a(A)', 'b(A)', 'c(A)', 'b/a', 'c/a', 'c_ext/a_ext']
    for i in range(nfigure):
        plt.subplot(nfigure, 1, i + 1)
        plt.ylabel(lp_ylabel[i])
        plt.plot(vol, lp_ori[i], 'o-', vv, lp_ext[i], 'rd')
        if i == 0:
            plt.legend(['ori', 'ext'], loc='best')
    plt.xlabel('Volume (A**3)')
    plt.savefig('ext_velp.png')
    if show_fig:
    	plt.show()
    plt.close()

    # save the fit data
    # get vba.dat
    fw = open(fout, 'w+')
    fw.write('#%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s\n' %
             ('volume', 'eng', 'cella', 'cellb', 'cellc', 'b/a', 'c/a', 'cext/aext'))
    for i in range(ndata):
        fw.write('%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\n' %
                 (vv[i], ee[i], cellaa[i], cellbb[i], cellcc[i], cellbaba[i], cellcaca[i], cellca_cal[i]))
        fw.flush()
    fw.close()
    print("\n>> Storing the extrapolate results in %s\n" % fout)
    print("\n>> DONE!")
    return


def lsqfit_eos(func, fin, par, fstart, fend, show_fig=False, fout='EoSfit.out', refit=-1):
    # make the screen output better.
    print('\n')
    print("\t>> We are using [ %s ] to fit the V-E relationship << \t" % func)
    print('\n')

    fs = fstart
    fe = fend
    p = par

    if fs < 0:
        print('start fitting range index must be a positive integer!')
        os.sys.exit(-1)

    # p0 = [e0, b0, bp, v0, bpp]
    if refit == -1:
        if (func == 'morse_AB') or (func == 'morse_3p'):
            A = 6
            B = 0.5 * A
            if func == 'morse_AB':
                p0 = [p[0], A, B, p[3]]
            else:
                p0 = [p[0], A, p[3]]
        elif func in ['mie', 'mie_simple']:
            p0 = [p[0], 4, 6, p[3]]
        elif func == 'morse_6p':
            P = 1
            Q = 2
            m = 1
            n = 1
            p0 = [p[0], p[1], P, p[3], Q, m, n]
        elif func == 'SJX_5p':
            alpha = 1
            beta = 1
            n = 1
            p0 = [p[0], alpha, beta, p[3], n]
        else:
            p0 = p
        print(
            '>> use initial guess of fitted-paramters by polynomial fitting results:\n')
    else:
        p0 = np.loadtxt(fout)
        print(
            '>> use initial guess of fitted-paramters by last fitting results:\n')
    print(p0)

    vol, en = read_ve(fin)
    ndata = len(vol)
    # To determine the data points included in the fitting
    if fe > ndata:
        fe = ndata
        print(
            '\n[WARNING]: f_end exceed available data numbers, Reset it to be %d now.' % ndata)
    if fs > ndata:
        print(
            '\n[WARNING]: f_start exceed available data numbers, Reset it to be 1: %d now.' % ndata)
        print('and Reset f_end to be %d now.' % ndata)
        fs = 1
        fe = ndata

    if fe == -1:
        print('\n[ATTENTIONS]: fend = -1')
        num = ndata - fs + 1
    else:
        num = fe - fs + 1

    vol = vol[fs - 1: fe]
    en = en[fs - 1: fe]
    print('\n%d/%d data was used in the fitting ...\n' % (num, ndata))

    #*************************************************************************
    # fit it: step 1.
    popt, pcov, infodict, mesg, ier = leastsq(
        eval('res_' + func), p0, args=(en, vol), full_output=1, maxfev=(len(p0) + 1) * 400)
    nfev = infodict['nfev']
    fvec = infodict['fvec']

    '''
    psi_i = sum(fvec**2)
    psi_min = min(fvec)
    pfactor = np.zeros(num)
    wi = (4 + 1) / num * psi_i / psi_min
    pfi = np.exp(-wi * wi)
    '''
    if ier not in [1, 2, 3, 4]:
        raise RuntimeError("Optimal parameters not found: " + mesg)
    #*************************************************************************

    print('*' * 80)
    print(">> fitted parameters (with %d iterations):" % nfev)
    # p0 = [e0, b0, bp, v0, bpp]
    if func == 'morse_AB':
        e0, A, B, v0 = popt
        print("%12s\t%12s\t%12s\t%12s" % ('V0(A**3)', 'A', 'B', 'E0(eV)'))
        print('%12f\t%12f\t%12f\t%12f\n' % (v0, A, B, e0))
    elif func == 'morse_3p':
        e0, A, v0 = popt
        B = 0.5 * A
        print("%12s\t%12s\t%12s\t%12s" % ('V0(A**3)', 'A', 'B', 'E0(eV)'))
        print('%12f\t%12f\t%12f\t%12f\n' % (v0, A, B, e0))
    elif func in ['mie', 'mie_simple']:
        e0, m, n, v0 = popt
        print("%12s\t%12s\t%12s\t%12s" % ('V0(A**3)', 'm', 'n', 'E0(eV)'))
        print('%12f\t%12f\t%12f\t%12f\n' % (v0, m, n, e0))
    elif func == 'morse_6p':
        e0, b0, bp, v0, bpp, m, n = calc_props_morse_6p(popt)
        b0 = eV2GPa * b0
        bpp = bpp / eV2GPa
        print("%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s" %
              ('V0(A**3)', 'B0(GPa)', 'Bp', 'E0(eV)', 'Bpp(1/GPa)', 'm', 'n'))
        print('%12f\t%12f\t%12f\t%12f\t%12f\t%12f\t%12f\n' %
              (v0, b0, bp, e0, bpp, m, n))
    elif func == 'SJX_5p':
        e0, b0, bp, v0, n = calc_props_SJX_5p(popt)
        b0 = eV2GPa * b0
        print("%12s\t%12s\t%12s\t%12s\t%12s" %
              ('V0(A**3)', 'B0(GPa)', 'Bp', 'E0(eV)', 'n'))
        print('%12f\t%12f\t%12f\t%12f\t%12f\n' % (v0, b0, bp, e0, n))
    elif func in ['mBM4poly', 'mBM5poly', 'mBM4', 'LOG4', 'vinet', 'morse', 'BM4']:
        prop_func = eval('calc_props_' + func)
        e0, b0, bp, v0, bpp = prop_func(popt)
        b0 = eV2GPa * b0
        bpp = bpp / eV2GPa
        print("%12s\t%12s\t%12s\t%12s\t%12s" %
              ('V0(A**3)', 'B0(GPa)', 'Bp', 'E0(eV)', 'Bpp(1/GPa)'))
        print('%12f\t%12f\t%12f\t%12f\t%12f\n' % (v0, b0, bp, e0, bpp))
    else:
        e0, b0, bp, v0, bpp = popt
        b0 = eV2GPa * b0
        bpp = bpp / eV2GPa
        print("%12s\t%12s\t%12s\t%12s\t%12s" %
              ('V0(A**3)', 'B0(GPa)', 'Bp', 'E0(eV)', 'Bpp(1/GPa)'))
        print('%12f\t%12f\t%12f\t%12f\t%12f\n' % (v0, b0, bp, e0, bpp))

    # write the fitted results in fit.out
    fw = open(fout, 'w+')
    for i in range(len(popt)):
        fw.write("%f\n" % popt[i])
        fw.flush()
    fw.close()

    # data evaluate by fitted parameters --to be used in the plot
    vol_i = np.linspace(min(vol) * 0.95, max(vol) * 1.05, len(vol) * 2)
    en_i = repro_ve(func, vol_i, popt)

    print('*' * 80)
    # calculate fitted residuals and square errors
    res_opt = np.zeros(len(vol))
    for i in range(len(vol)):
        res_opt[i] = (fvec[i])**2
    fit_res = sum(res_opt)
    fit_var = np.var(fvec)
    fit_std = np.std(fvec)
    print("\nfitted residuals\t= %16e\n" % fit_res)
    print("fitted variations\t= %16e\n" % fit_var)
    print("standard deviations\t= %16e\n" % fit_std)
    # if fit_res > 1e-4:
    #    print("\n>> Residuals seems too large, please refit it by swithing argument --refit 1!\n")
    #    show = 'F'  # reset show tag, not to show the figure.
    plt.plot(vol, en, 'o', vol_i, en_i)
    plt.title('EoS fitted by: %s model' % str(func))
    plt.legend(['calc', func + '-fit'], loc='best')
    plt.xlabel('Volume (A**3)')
    plt.ylabel('Energy (eV)')
    plt.savefig('EoSfit_' + func + '.png')
    if show_fig:
        plt.show()
    print('*' * 80)
    plt.close()

    # reproduce data by the popt and write it into files.
    repro_en = repro_ve(func, vol, popt)
    repro_press = repro_vp(func, vol, popt)
    fve = open(func + '_ve_fit.dat', 'w+')
    fvp = open(func + '_vp_fit.dat', 'w+')
    fve.write('#%20s\t%20s\t%20s\t%20s\n' %
              ('volume(A**3)', 'energy(fit)', 'energy(cal)', 'dE(%)'))
    fvp.write('#%20s\t%20s\t%20s\t%20s\n'
              % ('volume(A**3)', 'pressure(GPa)', 'pressure(Mbar)', 'pressure(kbar)'))
    for i in range(len(vol)):
        fve.write('%20f\t%20f\t%20f\t%20f\n' %
                  (vol[i], repro_en[i], en[i], 100 * np.abs((en[i] - repro_en[i]) / en[i])))
        fve.flush()
        p_tmp = repro_press[i]
        fvp.write('%20f\t%20f\t%20f\t%20f\n' %
                  (vol[i], p_tmp, p_tmp / 100, p_tmp * 10))
        fvp.flush()
    fve.close()
    fvp.close()
    return popt


def parse_argument():
    parser = argparse.ArgumentParser(
        description=" Script to fit EOS in MTP calculations")
    parser.add_argument('choice', help='to run mfp, ext_vec, ext_velp?')
    parser.add_argument(
        'infile', help='input ve.dat, vec.dat, velp.dat, or vlp.dat')
    parser.add_argument('-eos', '--eos', default='morse_AB',
                        help='Please choose one of the EOS name in: ' + str(get_eos_list()))
    parser.add_argument(
        '--show', default='T', help='to show the fitted results[T] or not[F]?')
    parser.add_argument('-vr', '--vrange', type=float, nargs=3,
                        default=[30, 50, 101], help="extend range used for ext, e.g. [vols, vole, ndata]: 1 10 31")
    parser.add_argument(
        '-fr', '--frange', type=int, nargs=2, default=[1, 1000], help='data range to be fitted, e.g. [ns, ne]')
    parser.add_argument('--refit', type=int, default=-1,
                        help='fit the data with last fitted parameters as initial guess')
    parser.add_argument('-eord', '--extorder', type=int,
                        default=3, help='spline order used in extrapolate')
    parser.add_argument('-enum', '--extnum', type=int, default=20,
                        help='number of data to be extrapolate')
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s ' + __version__())
    args = parser.parse_args()
    return(args)


# main
if __name__ == '__main__':
    args = parse_argument()
    cho = args.choice
    fin = args.infile
    func = args.eos
    SHOW = args.show
    REFIT = args.refit

    eos_list = get_eos_list()
    if func not in eos_list:
        print("Choose one of the following eos name:")
        print(eos_list)
        os.sys.exit(-1)
    else:
        pass
    # p0 = [e0, b0, bp, v0, bpp]

    # used in mfp choice
    fs = args.frange[0]
    fe = args.frange[1]

    if fs < 0:
        print("ERROR, start range index must be a positive integer!")
        os.sys.exit(-1)

    if cho == 'mfp':
        p0 = init_guess(fin)
        lsqfit_eos(
            func, fin, p0, fs, fe, SHOW, refit=REFIT)
    elif cho == 'ext_vec':
        # used in ext choice
        ExtVS = args.vrange[0]
        ExtVE = args.vrange[1]
        ExtNum = int(args.vrange[2])
        if ExtVS < 0 or ExtVE < 0 or ExtNum < 0:
            print("ERROR, range setting must be a positive value!")
            exit(1)
        # if ExtNum % 2 == 0:
        #    print("[WARNING]: ndata = %d, reset it to be %d" %
        #          (ExtNum, ExtNum + 1))
        #    ExtNum = ExtNum + 1
        p0 = init_guess(fin)
        ext_vec(func, fin, p0, fs, fe, ExtVS, ExtVE, ExtNum, refit=REFIT)
    elif cho == 'ext_velp':
        EORDER = args.extorder
        ExtVS = args.vrange[0]
        ExtVE = args.vrange[1]
        ExtNum = int(args.vrange[2])
        if ExtVS < 0 or ExtVE < 0 or ExtNum < 0:
            print("ERROR, range setting must be a positive value!")
            exit(1)
        if ExtNum % 2 == 0:
            print("[WARNING]: ndata = %d, reset it to be %d" %
                  (ExtNum, ExtNum + 1))
            ExtNum = ExtNum + 1
        ext_velp(fin, fs, fe, ExtVS, ExtVE, ExtNum,
                 order=EORDER, method='unispl')
    else:
        print('Unkown Choice, abort now ... ...')
        os.sys.exit(-1)
