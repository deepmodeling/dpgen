#----Force field------------------------------------------------
neigh_modify	every 1 delay 0 check no
pair_style	deepmd graph.000.pb graph.001.pb graph.002.pb graph.003.pb 10 model_devi.out
#pair_style	deepmd frozen_model.pb
pair_coeff


