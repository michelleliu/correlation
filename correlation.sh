skip=0
total_time=20000
meas_time=1000
sep=13.0
inp_file=../../npt_13.0/C-walls_13.0_npt_1atm_298K.t28502.traj
suffix=confined
#./correlation -s $skip -t $total_time -m $meas_time -d $sep -e $suffix $inp_file

inp_file=../../clean_npt_water/neat_water.t28505.traj
suffix=neat
./correlation -s $skip -t $total_time -m $meas_time -d $sep -e $suffix $inp_file
