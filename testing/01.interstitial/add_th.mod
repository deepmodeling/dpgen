#----Add interstitial atom--------------------------------------
variable	mhao equal ${ao}*1.0/4.0 
create_atoms	1 single ${mhao} ${mhao} ${mhao} units box
variable	r2 equal ${ao}/100.0
region		select sphere ${mhao} ${mhao} ${mhao} ${r2} units box
group		fnull region select
