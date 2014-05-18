skip=0
total_time=20000
meas_time=1000
sep=13.0
inp_file=../../npt_13.0/C-walls_13.0_npt_1atm_298K.t28502.traj
O_id=2

#------------
#suffix=full

#inp_file=../../clean_npt_water/neat_water.t28505.traj
#suffix=neat

#inp_file=../npt_water/water_sandia.t28615.traj
#suffix=sandia

inp_file=../npt_water/orsi/TIP4P-Ew.t28634.traj
suffix=tip4pew
./correlation-2 -s $skip -t $total_time -m $meas_time -d $sep -e $suffix $inp_file

# inp_file=../npt_water/orsi/TIP3P-Ew.t28633.traj
# suffix=tip3pew
#
# total_time=5000
# meas_time=1000
# inp_file=../npt_water/MICHELLE.txt-cleanup
# suffix=spce-tinker

#./correlation-tinker -s $skip -t $total_time -m $meas_time -d $sep -e $suffix $inp_file
