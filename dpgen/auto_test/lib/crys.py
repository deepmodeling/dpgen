import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure


def fcc(ele_name="ele", a=4.05):
    box = np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
    box *= a
    return Structure(box, [ele_name], [[0, 0, 0]])


def fcc1(ele_name="ele", a=4.05):
    latt = Lattice.cubic(a)
    return Structure(
        latt,
        [ele_name, ele_name, ele_name, ele_name],
        [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
    )


def sc(ele_name="ele", a=2.551340126037118):
    latt = Lattice.cubic(a)
    return Structure(latt, [ele_name], [[0, 0, 0]])


def bcc(ele_name="ele", a=3.2144871302356037):
    latt = Lattice.cubic(a)
    return Structure(
        latt,
        [ele_name, ele_name],
        [
            [0, 0, 0],
            [0.5, 0.5, 0.5],
        ],
    )


def hcp(
    ele_name="ele", a=4.05 / np.sqrt(2), c=4.05 / np.sqrt(2) * 2.0 * np.sqrt(2.0 / 3.0)
):
    box = np.array([[1, 0, 0], [0.5, 0.5 * np.sqrt(3), 0], [0, 0, 1]])
    box[0] *= a
    box[1] *= a
    box[2] *= c
    latt = Lattice(box)
    return Structure(
        latt, [ele_name, ele_name], [[0, 0, 0], [1.0 / 3, 1.0 / 3, 1.0 / 2]]
    )


def dhcp(
    ele_name="ele", a=4.05 / np.sqrt(2), c=4.05 / np.sqrt(2) * 4.0 * np.sqrt(2.0 / 3.0)
):
    box = np.array([[1, 0, 0], [0.5, 0.5 * np.sqrt(3), 0], [0, 0, 1]])
    box[0] *= a
    box[1] *= a
    box[2] *= c
    latt = Lattice(box)
    return Structure(
        latt,
        [ele_name, ele_name, ele_name, ele_name],
        [
            [0, 0, 0],
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 4.0],
            [0, 0, 1.0 / 2.0],
            [2.0 / 3.0, 2.0 / 3.0, 3.0 / 4.0],
        ],
    )


def diamond(ele_name="ele", a=2.551340126037118):
    box = np.array([[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]])
    box *= a
    latt = Lattice(box)
    return Structure(
        latt,
        [ele_name, ele_name],
        [
            [0.12500000000000, 0.12500000000000, 0.12500000000000],
            [0.87500000000000, 0.87500000000000, 0.87500000000000],
        ],
    )
