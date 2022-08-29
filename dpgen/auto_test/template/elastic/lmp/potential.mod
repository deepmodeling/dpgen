# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# ================= Choose potential ========================
pair_style	deepmd graph.000.pb graph.001.pb graph.002.pb graph.003.pb  400 model_devi.out
#pair_style	deepmd frozen_model.pb
pair_coeff     * *

# Setup neighbor style
neigh_modify	every 1 delay 0 check yes

# Setup minimization style
min_style	cg
min_modify	dmax ${dmax} 

# Setup output
thermo		10
thermo_style	custom step ke pe press pxx pyy pzz pxy pxz pyz lx ly lz
thermo_modify norm no
