#----Add interstitial atom--------------------------------------
# r2 is the radius of the copper atom
# variable	r2 equal sqrt(${ao}^2+${ao}^2)/4
# region		select sphere 0 0 0 ${r2} units box
variable	r2 equal ${ao}/10.0
region		del sphere 0 0 0 ${r2} units box
delete_atoms	region del
variable	qhao_0 equal ${ao}/6.0 
variable	qhao_1 equal -${ao}/6.0 
create_atoms	1 single ${qhao_0} ${qhao_0} 0 units box
create_atoms	1 single ${qhao_1} ${qhao_1} 0 units box

