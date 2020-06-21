import os, glob, warnings
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
import numpy as np
import dpdata


def make_repro(vasp_lmp_path, path_to_work):
    if not os.path.exists(vasp_lmp_path):
        raise RuntimeError("please do VASP or LAMMPS calcualtions first")
    vasp_lmp_task = glob.glob(os.path.join(vasp_lmp_path, 'task.[0-9]*[0-9]'))
    assert len(vasp_lmp_task) > 0, "Please do VASP or LAMMPS calculations first"
    vasp_lmp_task.sort()
    task_num = 0
    task_list = []
    for ii in vasp_lmp_task:
        # get vasp or lmp energy
        outcar = os.path.join(ii, 'OUTCAR')
        log_lmp = os.path.join(ii, 'log.lammps')
        if os.path.exists(outcar):
            energies = vasp.get_energies(outcar)
            # get xdat
            xdatcar = os.path.join(ii, 'XDATCAR')
            os.chdir(path_to_work)
            if os.path.exists('XDATCAR'):
                os.remove('XDATCAR')
            os.symlink(os.path.relpath(xdatcar), 'XDATCAR')
            xdat_lines = open('XDATCAR', 'r').read().split('\n')
            natoms = vasp.poscar_natoms('XDATCAR')
            xdat_secsize = natoms + 8
            xdat_nframes = len(xdat_lines) // xdat_secsize
            if xdat_nframes > len(energies):
                warnings.warn('nframes %d in xdatcar is larger than energy %d, use the last %d frames' % (
                    xdat_nframes, len(energies), len(energies)))
                xdat_nlines = -1 * len(energies) * xdat_secsize  # 06/12 revised
                xdat_lines = xdat_lines[xdat_nlines:]
            xdat_nframes = len(xdat_lines) // xdat_secsize
            print(xdat_nframes, len(energies))

            # loop over frames
            for jj in range(xdat_nframes):
                output_task = os.path.join(path_to_work, 'task.%06d' % task_num)
                task_num += 1
                task_list.append(output_task)
                os.makedirs(output_task, exist_ok=True)
                os.chdir(output_task)
                # clear dir
                for kk in ['INCAR', 'POTCAR', 'POSCAR.orig', 'POSCAR', 'conf.lmp', 'in.lammps']:
                    if os.path.exists(kk):
                        os.remove(kk)
                # make conf
                with open('POSCAR', 'w') as fp:
                    fp.write('\n'.join(xdat_lines[jj * xdat_secsize:(jj + 1) * xdat_secsize]))

        elif os.path.exists(log_lmp):
            # get dump
            dump = os.path.join(ii, 'dump.relax')
            if not os.path.exists(dump):
                raise FileNotFoundError("the LAMMPS calculations should output the snapshots as dump.relax")

            energies = get_energy_fordump(dump, log_lmp)

            # loop over frames
            for jj in range(len(energies)):
                output_task = os.path.join(path_to_work, 'task.%06d' % task_num)
                task_num += 1
                task_list.append(output_task)
                os.makedirs(output_task, exist_ok=True)
                os.chdir(output_task)
                # clear dir
                for kk in ['INCAR', 'POTCAR', 'POSCAR.orig', 'POSCAR', 'conf.lmp', 'in.lammps']:
                    if os.path.exists(kk):
                        os.remove(kk)

                # make conf
                d_dump = dpdata.System(dump, fmt='lammps/dump')
                d_dump.to('vasp/poscar', 'POSCAR', frame_idx=jj)

    return task_list


def get_energy_fordump(dump, log_lmp):
    dump = os.path.basename(dump)
    with open(log_lmp, 'r') as fin:
        lines = fin.read().split('\n')
    for idx, jj in lines:
        if dump == jj.split()[5]:
            dump_freq = int(jj.split()[4])
        if 'Step PotEng' in jj:
            step_idx = idx
        if 'Loop time' in jj:
            end_idx = idx
    energies = []

    for jj in lines[step_idx + 1, end_idx]:
        step = int(jj.split()[0])
        if step % dump_freq == 0:
            energies.append(float(jj.split()[1]))

    return energies


def post_repro(vasp_lmp_path, all_tasks, ptr_data):
    ptr_data += "Reproduce: Initial_path DFT_E(eV/atom)  LMP_E(eV/atom)  Difference(eV/atom)\n"
    vasp_lmp_task = glob.glob(os.path.join(vasp_lmp_path, 'task.[0-9]*[0-9]'))
    assert len(vasp_lmp_task) > 0, "Please do VASP or LAMMPS calcualtions first"
    vasp_lmp_task.sort()
    vasp_ener_tot = []
    lmp_ener_tot = []
    res_data = {}

    for ii in vasp_lmp_task:
        # compute vasp or lammps
        outcar = os.path.join(ii, 'OUTCAR')
        loglammps = os.path.join(ii, 'log.lammps')
        if os.path.exists(outcar):
            vasp_ener = np.array(vasp.get_energies(outcar))
            vasp_ener_file = os.path.join(ii, 'ener.vasp.out')
            # compute reprod
            lmp_ener = []

            if len(all_tasks) < (len(vasp_ener_tot) + len(vasp_ener)):
                raise RuntimeError("lammps tasks reproduced not equal to vasp")

            natoms = 1
            for jj in range(len(vasp_ener_tot), (len(vasp_ener_tot) + len(
                    vasp_ener))):  # all_tasks[len(vasp_ener_tot):(len(vasp_ener_tot) + len(vasp_ener))]:
                log_lmp = os.path.join(all_tasks[jj], 'log.lammps')
                if not os.path.exists(log_lmp):
                    raise RuntimeError("lammps reproduce not finished")
                natoms, epa, vpa = lammps.get_nev(log_lmp)
                lmp_ener.append(epa)
                lmp_ener_tot.append(epa)
                vasp_epa = list(vasp_ener)[jj - len(vasp_ener_tot)] / natoms
                ptr_data += '%s %7.3f  %7.3f  %7.3f\n' % (ii, vasp_epa,
                                                          epa, epa - vasp_epa)
            lmp_ener = np.array(lmp_ener)
            lmp_ener = np.reshape(lmp_ener, [-1, 1])
            vasp_ener_tot += list(vasp_ener)
            vasp_ener = np.reshape(vasp_ener, [-1, 1]) / natoms
            error_start = 1
            lmp_ener -= lmp_ener[-1] - vasp_ener[-1]
            diff = lmp_ener - vasp_ener
            diff = diff[error_start:]
            error = np.linalg.norm(diff) / np.sqrt(np.size(lmp_ener)-1)   # start from diff[1], so one element fewer
            res_data[ii] = {'nframes': len(vasp_ener), 'error': error}
            np.savetxt(vasp_ener_file, vasp_ener[error_start:])

        elif os.path.exists(loglammps):
            dump = os.path.join(ii, 'dump.relax')
            lmp_ener = np.array(get_energy_fordump(dump, loglammps))
            lmp_ener_file = os.path.join(ii, 'ener.lammps.out')
            # compute reprod
            vasp_ener = []

            if len(all_tasks) < (len(lmp_ener_tot) + len(lmp_ener)):
                raise RuntimeError("vasp tasks reproduced not equal to lammps")

            for jj in range(len(lmp_ener_tot), (len(lmp_ener_tot) + len(lmp_ener))):
                outcar = os.path.join(all_tasks[jj], 'OUTCAR')
                if not os.path.exists(outcar):
                    raise RuntimeError("vasp reproduce not finished")
                natoms, epa, vpa = vasp.get_nev(outcar)
                vasp_ener.append(epa)
                vasp_ener_tot.append(epa)
                lmp_epa = list(lmp_ener)[jj - len(lmp_ener_tot)] / natoms
                ptr_data += '%s %7.3f  %7.3f  %7.3f\n' % (ii, epa,
                                                          lmp_epa, lmp_epa - epa)

            vasp_ener = np.array(vasp_ener)
            vasp_ener = np.reshape(vasp_ener, [-1, 1])
            lmp_ener_tot += list(lmp_ener)
            lmp_ener = np.reshape(lmp_ener, [-1, 1]) / natoms
            # counting error from the first frame
            error_start = 0
            lmp_ener -= lmp_ener[-1] - vasp_ener[-1]
            diff = lmp_ener - vasp_ener
            diff = diff[error_start:]
            error = np.linalg.norm(diff) / np.sqrt(np.size(lmp_ener))
            res_data[ii] = {'nframes': len(lmp_ener), 'error': error}
            np.savetxt(lmp_ener_file, lmp_ener[error_start:])

    if not len(vasp_ener_tot) == len(lmp_ener_tot):
        raise RuntimeError("lammps tasks reproduced not equal to vasp")
    #    for ii in range(len(lmp_ener_tot)):
    #        ptr_data += '%7.3f  %7.3f  %7.3f\n' % (vasp_ener_tot[ii], lmp_ener_tot[ii],
    #                                               lmp_ener_tot[ii] - vasp_ener_tot[ii])
    return res_data, ptr_data
