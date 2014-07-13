skip=0
total_time=$1 #4000000
meas_time=4000 #50000
num_particles=1728

sep=13.0
#inp_file=../npt_13.0/q_C-walls_13.0_npt_1atm_298K_nnm.t32674.traj
inp_file=../npt_13.0/q-compute_58731.out
inp_file=../npt_13.0/mu-13.0-4ns-1728wat.out
suffix=c13-4ns
#./bin/correlation -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles $inp_file
#./bin/dipole -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles $inp_file

#inp_file=../npt_9.8/q_C-walls_9.8-npt_1atm_298K.t36084.traj
inp_file=../npt_9.8/q-compute_9.8_58737
inp_file=../npt_9.8/mu-9.8-4ns-1728wat.out
suffix=c9.8-4ns
sep=9.8

#inp_file=../npt_6.8/q_C-walls_6.8_fix-rigid-npt_1atm_298K.t36091.traj
inp_file=../npt_6.8/mu-6.8-4ns-1728wat.out
suffix=c6.8-4ns
sep=6.8

#inp_file=../npt_water/q_compute_36185.log-cleanup
inp_file=../npt_water/hbonds_1728-SPCE.t106473.traj
#inp_file=../npt_water/min_216-SPCE.t106451.traj
suffix=spce-4ns
num_particles=1728
sep=13
oid=2
./bin/hbonds -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles -o $oid -C $inp_file
