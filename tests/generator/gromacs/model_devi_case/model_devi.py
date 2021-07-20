#!/data1/anguse/zijian/deepmd-kit-devel/virtual_tf_2/bin/python
import deepmd.DeepPot as DP
import dpdata
import numpy as np
import os, sys
import json
def calc_model_devi_f(fs):
    '''
        fs : numpy.ndarray, size of `n_models x n_frames x n_atoms x 3`
    '''
    fs_mean = np.mean(fs, axis=0)
    # print(fs_mean.shape)
    fs_err = np.sum((fs - fs_mean) ** 2, axis=-1)
    # print(fs_err.shape)
    fs_devi = np.mean(fs_err, axis=0) ** 0.5
    # print(fs_devi.shape) 
    max_devi_f = np.max(fs_devi, axis=1)
    # min_devi_f = np.min(fs_devi, axis=1)
    # avg_devi_f = np.mean(fs_devi, axis=1)
    return max_devi_f

def write_model_devi_out(system, models, fname=None, trj_freq = 10):
    forces = []
    for model in models:
        labeled = system.predict(model)
        forces.append(labeled['forces'])
    forces = np.array(forces)
    max_devi_f = calc_model_devi_f(forces)
    model_devi_out = np.zeros((system.get_nframes(), 7))
    model_devi_out[:, 0] += np.arange(system.get_nframes()) * trj_freq
    model_devi_out[:, 4] += max_devi_f
    if fname is not None:
        np.savetxt(fname,
                   model_devi_out,
                   fmt=['%d'] + ['%.8e' for _ in range(6)],
                   delimiter='\t',
                   header='step\tmax_devi_e\tmin_devi_e\tavg_devi_e\tmax_devi_f\tmin_devi_f\tavg_devi_f')
    return model_devi_out

if __name__ == "__main__":
    system = dpdata.System(sys.argv[1], fmt='gromacs/gro')
    if os.path.isfile("job.json"):
        trj_freq = json.load(open("job.json")).get("trj_freq", 10)
    else:
        trj_freq = 10
    if not os.path.isdir("traj"):
        os.mkdir("traj")
    for i in range(system.get_nframes()):
        system[i].to_gromacs_gro("traj/%d.gromacstrj" % (trj_freq * i) )
    models = [DP(f"../graph.{ii:03}.pb") for ii in range(4)]
    write_model_devi_out(system, models, "model_devi.out", trj_freq)
