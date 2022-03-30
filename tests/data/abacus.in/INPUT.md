INPUT_PARAMETERS

calculation md
atom_file           STRU
kpoint_file         KPT
pseudo_dir ./
ntype 2
symmetry 1
ecutwfc 90

npool 4
charge_extrap second-order
mixing_type pulay
mixing_beta 0.4
dr2 1e-6
nstep 1
force_thr_ev 0.02
move_method cg
out_stru 0
