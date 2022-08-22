INPUT_PARAMETERS

calculation md
pseudo_dir ./
ntype 2
symmetry 1
ecutwfc 90
scf_thr 1e-2
mixing_type pulay
mixing_beta 0.4
basis_type lcao
md_nstep 3
cal_stress 1
deepks_model model.ptg

md_tfirst 10
md_tfreq 0.1
