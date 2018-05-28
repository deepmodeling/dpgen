#----Add interstitial atom--------------------------------------
# r2 is the radius of the copper atom
# variable	r2 equal sqrt(${ao}^2+${ao}^2)/4
# region		select sphere 0 0 0 ${r2} units box
variable	mhao equal -${ao}/2.0 
create_atoms	1 single 0 ${mhao} 0 units box
