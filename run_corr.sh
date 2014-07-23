skip=0
total_time=$1 #4000000
meas_time=$2 #50000
num_particles=1728

inp_file=../npt_13.0/mu_xyzq_c13.0.t142348.traj
sep=13.0
suffix=c13
oid=2
#./bin/hbonds -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles -o $oid -C $inp_file
#./bin/correlation -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles -o $oid $inp_file
inp_file=../npt_13.0/mu_xyzq_c13.0.m142348.out
#./bin/dipole -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles $inp_file

inp_file=../npt_9.8/mu_xyzq_c9.8.t142345.traj
suffix=c9.8
oid=2
sep=9.8
./bin/hbonds -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles -o $oid -C $inp_file
#./bin/correlation -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles -o $oid $inp_file
inp_file=../npt_9.8/mu_xyzq_c9.8.m142345.out
#./bin/dipole -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles $inp_file

inp_file=../npt_6.8/mu_xyzq_c6.8.t142344.traj
suffix=c6.8
oid=2
sep=6.8
#./bin/hbonds -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles -o $oid -C $inp_file
#./bin/correlation -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles -o $oid $inp_file
inp_file=../npt_6.8/mu_xyzq_c6.8.m142344.out
#./bin/dipole -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles $inp_file

#inp_file=../npt_water/q_compute_36185.log-cleanup
inp_file=../npt_water/hbonds_1728-SPCE.t106475.traj
#inp_file=../npt_water/min_216-SPCE.t106451.traj
suffix=spce-4ns
num_particles=1728
sep=13
oid=2


inp_file=../coarse_13.0/coarse_walls_13.0_npt_1atm_298K.t142167.traj
num_particles=288
sep=13.0
oid=1
suffix=spce-coarse-13.0
#./bin/correlation -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles -o $oid $inp_file

inp_file=../coarse_9.8/coarse_walls_9.8_npt_1atm_298K.t142165.traj
suffix=spce-coarse-9.8
num_particles=288
sep=9.8
oid=1
#./bin/correlation -s $skip -t $total_time -m $meas_time -d $sep -e $suffix -w $num_particles -o $oid $inp_file
