INPUT_PARAMETERS
#Parameters (1.General)
suffix			ABACUS
calculation             md
ntype			1
symmetry		0

#Parameters (2.Iteration)
ecutwfc			100
scf_thr			1e-8
scf_nmax		100
#Parameters (3.Basis)
basis_type		lcao

#Parameters (4.Smearing)
smearing_method		gauss
smearing_sigma			0.002

#Parameters (5.Mixing)
mixing_type		pulay
mixing_beta		0.3

cal_stress 1

md_type 	0
md_nstep        10
md_tfirst 10
md_tfreq        0.5
