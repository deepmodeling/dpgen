INPUT_PARAMETERS
#Parameters (1.General)
suffix			ABACUS
calculation     md
ntype			1
symmetry		1
pseudo_dir      ./
pseudo_type     upf201

#Parameters (2.Iteration)
ecutwfc			40
scf_thr				1e-6
scf_nmax			100


#Parameters (3.Basis)
basis_type		pw

#Parameters (4.Smearing)
smearing_method		gauss
smearing_sigma			0.002

#Parameters (5.Mixing)
mixing_type	    broyden 
mixing_beta		0.7

cal_stress 1
md_nstep 3
md_tfirst 10
md_tfreq 0.1
out_stru 1
