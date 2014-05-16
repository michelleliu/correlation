skip=0
total_time=20000
meas_time=1000
sep=13.0
#inp_file=../npt_6.8/C-walls_6.8_fix-rigid-npt_1atm_298K.t11136.traj
inp_file=../C-walls_13.0_npt_1atm_298K.t28502.traj

#python correlation.py -s $skip -t $total_time -m $meas_time -d $sep $inp_file
./correlation -s $skip -t $total_time -m $meas_time -d $sep $inp_file

skip=0
total_time=20000
meas_time=1000
sep=13.0
#inp_file=../npt_13.0/C-walls_13.0_fix-rigid-npt_1atm_298K.t11169.traj

#python correlation-window.py -s $skip -t $total_time -m $meas_time -d $sep $inp_file
#python correlation.py -s $skip -t $total_time -m $meas_time -d $sep $inp_file
