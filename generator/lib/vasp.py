#!/usr/bin/python3 

def make_vasp_kpoints (kpoints) :
    ret = ""
    ret += "Automatic mesh\n"
    ret += "0\n"
    ret += "Gamma\n"
    ret += "%d %d %d\n" % (kpoints[0], kpoints[1], kpoints[2])
    ret += "0  0  0\n"
    return ret

def make_vasp_incar (ecut, ediff, npar, kpar, kspacing = 0.5, kgamma = True) :
    ret = ''
    ret += 'PREC=A\n'
    ret += 'ENCUT=%d\n' % ecut
    ret += 'ISYM=0\n'
    ret += 'ALGO=fast\n'
    ret += 'EDIFF=%e\n' % ediff
    ret += 'LREAL=T\n'
    ret += 'NPAR=%d\n' % npar
    ret += 'KPAR=%d\n' % kpar
    ret += "\n"
    ret += 'NELMIN=4\n'
    ret += 'ISIF=2\n'
    ret += 'ISMEAR=1\n'
    ret += 'SIGMA=0.25\n'
    ret += 'IBRION=-1\n'
    ret += "\n"
    ret += 'NSW=0\n'
    ret += "\n"
    ret += 'LWAVE=F\n'
    ret += 'LCHARG=F\n'
    ret += 'PSTRESS=0\n'
    ret += "\n"
    ret += 'KSPACING=%f\n' % kspacing
    if kgamma:
        ret += 'KGAMMA=.TRUE.\n'
    else :
        ret += 'KGAMMA=.FALSE.\n'

    return ret

def make_vasp_kpoints_gamma (kpoints) :
    ret = ''
    ret += 'Automatic mesh\n'
    ret += '0\n'
    ret += 'Gamma\n'
    ret += '%d %d %d\n' % (kpoints[0], kpoints[1], kpoints[2])
    ret += '0  0  0\n'
    return ret

def make_vasp_kpoints (kpoints) :
    return make_vasp_kpoints_gamma(kpoints)
