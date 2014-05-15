import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import math
import time
# example usage:
# python correlation.py -s 5000 -t 3000 -m 2000 -d 6.8 npt_6.8/rst63000.C-walls_6.8_fix-rigid-npt_rst.t10814.traj

def is_in_box(x,y,z,x0,x1,y0,y1,z0,z1):
    return x > x0 and x < x1 and y > y0 and y < y1 and z > z0 and z < z1

def sq_displacement(a,b,c,x,y,z):
    return (x-a)**2+(y-b)**2+(z-c)**2

def read_input_traj():
    N_occu_total = 0
    theta = [[0 for x in xrange(num_particles)] for x in xrange(num_timesteps)]
    if VELOCITY:
        omega_total = 0
        velocities = [[[0 for x in xrange(3)] for x in xrange(num_particles)] for x in xrange(num_timesteps)]
    if MSD:
        positions = [[[0 for x in xrange(3)] for x in xrange(num_particles)] for x in xrange(num_timesteps)]

    # skip a number of steps from beginning
    data.seek(0)
    for i in np.arange(skip_steps):
        for i in np.arange(num_atoms+9):
            data.readline()

    # beginning data parsing
    print "{0} Atoms, {1} Waters".format(num_atoms,num_particles)
    for step in np.arange(num_timesteps):
        if step % (num_timesteps/10) == 0:
            percent = 100.0*step/num_timesteps
            print "...{0}% done...".format(percent)
        #--- skipping non-xyz lines
        data.readline()
        current_time=data.readline()
        data.readline()
        data.readline()
        data.readline()
        data.readline() # box_x_bounds =
        data.readline() # box_y_bounds =
        data.readline() # box_z_bounds =
        data.readline()
        #---
        #num_waters_check = 0 # for testing
        # for int i in np.arange(2):
        #     box_x_bounds[i] = float(box_x_bounds[i])
        #     box_y_bounds[i] = float(box_y_bounds[i])
        #     box_z_bounds[i] = float(box_z_bounds[i])
        if PRINTING:
            print "Time: {0}".format(current_time)
            print "Computing step {0}...".format(step)
        #t = step*dump_interval + starttime
        particle_idx = 0
        for i in np.arange(num_atoms):
            line = data.readline()
            s = line.split()
            if float(s[2])==2: # !!! 0 if test.data or xyz, 2 if traj
                #num_waters_check += 1
                particle_x = float(s[3]) # !!! 1 if xyz, 2 if test.data, 3 if traj
                particle_y = float(s[4]) # !!!
                particle_z = float(s[5]) # !!!
                particle_idx = int(s[1])-1 # !!! only for traj ?
                particle_ix = int(s[9])
                particle_iy = int(s[10])
                particle_iz = int(s[11])
                if VELOCITY:
                    velocities[step][particle_idx][0] = float(s[6]) # vx
                    velocities[step][particle_idx][1] = float(s[7]) # vy
                    velocities[step][particle_idx][2] = float(s[8]) # vz
                    #print "v_x for particle {0} is {1}".format(particle_idx,velocities[step][particle_idx][0])
                # check if particle is in the intersolute space
                if MSD:
                    positions[step][particle_idx][0] = particle_x # x-coord
                    positions[step][particle_idx][1] = particle_y # y-coord
                    positions[step][particle_idx][2] = particle_z # z-coord
                if is_in_box(particle_x,particle_y,particle_z,x0,x1,y0,y1,z0,z1):
                    #print "I'm in the box! My coords: {0} {1} {2}, my ID: {3}".format(particle_x,particle_y, particle_z, particle_idx) # for testing
                    # if particle is in box, 1 the entry corresponding to that particle
                    theta[step][particle_idx] = 1
                    # add one to running sum of number of particles in space
                    N_occu_total+=1
                    # increment omega_total
                    if VELOCITY:
                        for k in np.arange(3):
                            omega_total += (velocities[step][particle_idx][k]*velocities[step][particle_idx][k])

        # this loop is where we calculate all the correlation functions
        if PRINTING:
            print "measuring correlation at step {0}".format(step)
        for step_meas in np.arange(min(step,len(R))):
            occu_tmp = 0
            vel_tmp = 0
            sq_disp_tmp = 0

            # initial time
            t0 = step-step_meas

            for j in np.arange(num_particles):
                # don't do anything if the particle wasn't in the box
                # at time t0
                if theta[t0][j]:

                    # occupation time distribution:
                    # find number of the particles in box at t=step
                    # who were also in the box step_meas steps ago
                    occu_tmp += theta[step][j]

                    # velocity autocorrelation
                    if VELOCITY:
                        # calculate dot product of velocities
                        # and add this dot product to the summand
                        for k in np.arange(3):
                            vel_tmp += (velocities[step][j][k]*velocities[t0][j][k])
                    if MSD:
                        # calculate the sq displacement
                        # and add this to the sum
                        sq_disp_tmp += ( sq_displacement(positions[step][j][0],positions[step][j][1],
                            positions[step][j][2],positions[t0][j][0],
                            positions[t0][j][1],positions[t0][j][2]) )


            if PRINTING:
                print "Found {0} correlations from {1} steps ago".format(occu_tmp, step_meas)
            # increments the summand
            R[step_meas] += occu_tmp
            if VELOCITY:
                Cv[step_meas] += vel_tmp
            if MSD:
                msdisp[step_meas] += sq_disp_tmp

        #print "{0} waters".format(num_waters_check)
        if PRINTING:
            print "----------------\n"

    N_occu_ave = N_occu_total*1.0/num_timesteps
    print "Pre-normalization R: {0}".format(R[:30])
    if PRINTING:
        print "N_occu_total = {0}, N_occu_ave = {1}".format(N_occu_total,N_occu_ave)
    if N_occu_ave == 0:
        print "N_occu_ave = 0; will not continue normalization"
    else:
        for i in np.arange(len(R)):
            t0_max = num_timesteps-i
            # normalize R
            R[i] = R[i]*1.0/N_occu_ave/(t0_max)
            if VELOCITY:
                # normalize Cv
                omega = 1.0*omega_total/N_occu_ave/(num_timesteps) # TODO: check if right t0max
                Cv[i] = Cv[i]*1.0/N_occu_ave/(t0_max)/omega
            if MSD:
                # normalize MSD
                msdisp[i] = msdisp[i]/(t0_max)

    if VELOCITY:
        print "N_occu_ave: {0}\nR: {1}\nCv: {2}".format(N_occu_ave,R[:30],Cv[:30])
    else:
        print "N_occu_ave: {0}\nR: {1}\n".format(N_occu_ave,R[:30])

def plot_rc(figname):
    fig1 = plt.figure()
    plt.plot(plot_time,R)
    plt.title("Occupation Time Distribution")
    plt.xlabel("time (ps)")
    plt.ylabel("$C_{R}(t)$")
    plt.draw()
    plt.savefig(figname)

def plot_cv(figname):
    fig2 = plt.figure()
    plt.plot(plot_time,Cv)
    plt.title("Velocity Autocorrelation Function")
    plt.xlabel("time (ps)")
    plt.ylabel("$C_{V}(t)$")
    plt.draw()
    plt.savefig(figname)

inputfile = sys.argv[len(sys.argv)-1]

data = open(inputfile)
line = data.readline()
starttime=int(data.readline())
data.readline()
num_atoms = int(data.readline())
PRINTING=False
VELOCITY=True
MSD=True
num_particles=1728

# set dimensions of observation box
x0=-6.05+15
x1=6.05+15
y0=-5.25+15
y1=5.25+15
z0=16.1
z1=22.90
# x0=-6.05
# x1=16.05
# y0=-5.25
# y1=15.25
# z0=15.8
# z1=43.25
num_timesteps=4000 # !!! length of time to collect data
if '-d' in sys.argv:
    z_center = 19.5
    z_dist = float(sys.argv[sys.argv.index('-d')+1])
    z0 = z_center - z_dist/2
    z1 = z_center + z_dist/2
z_dist = z1-z0
if '-t' in sys.argv:
    num_timesteps = int(sys.argv[sys.argv.index('-t')+1])
delta_t = 2.0 # fs
dump_interval = 10.0 # !!!
dump_time = dump_interval*delta_t

skip_steps = 10 # discard first 10 steps by default
if '-s' in sys.argv:
    skip_steps = int(sys.argv[sys.argv.index('-s')+1])

max_meas_time=1000 # this is how long the X-axis of the plot will be by default
                   # measured in dump steps
                   # ideally: [max_meas_time] = [desired_time]*[conversion]/[dump_interval]
if '-m' in sys.argv:
    max_meas_time = int(sys.argv[sys.argv.index('-m')+1])
RCoutputfile = "R_C_d{3}_s{0}_m{1}_t{2}_ocorr".format(skip_steps,max_meas_time,num_timesteps,z_dist)
CVoutputfile = "C_V_d{3}_s{0}_m{1}_t{2}_ocorr".format(skip_steps,max_meas_time,num_timesteps,z_dist)
logoutputfile = "log_d{3}_s{0}_m{1}_t{2}_ocorr".format(skip_steps,max_meas_time,num_timesteps,z_dist)
MSDoutputfile = "MSD_d{3}_s{0}_m{1}_t{2}".format(skip_steps,max_meas_time,num_timesteps,z_dist)
if '-w' in sys.argv:
    RCoutputfile = sys.argv[sys.argv.index('-w')+1]
RCout = open(RCoutputfile, "w+")
if VELOCITY:
    CVout = open(CVoutputfile, "w+")
if MSD:
    MSDout = open(MSDoutputfile, "w+")
logout = open(logoutputfile, "w+")
print "Discarding first {0} steps".format(skip_steps)
logout.write("Discarding first {0} steps\n".format(skip_steps))
print "Reading {0} steps".format(num_timesteps)
logout.write("Reading {0} steps\n".format(num_timesteps))
print "Measuring correlations up to {0} steps".format(max_meas_time)
logout.write("Measuring correlations up to {0} steps\n".format(max_meas_time))
print "reading " + inputfile + "...\n"
logout.write("reading {0}...\n".format(inputfile))
logout.write("Observation box: x ({0}, {1}), y ({2}, {3}), z ({4}, {5}); d = {6} Angstroms\n".format(x0,x1,y0,y1,z0,z1,z_dist))
R = np.zeros(min(num_timesteps,max_meas_time))
if VELOCITY:
    Cv = np.zeros(min(num_timesteps,max_meas_time))
if MSD:
    msdisp = np.zeros(max_meas_time+1)

ptime = time.time()
read_input_traj()
ptime = time.time() - ptime
print "time to analyze input traj: {0} s".format(ptime)

plot_time = np.zeros(len(R))
# fix time axes
for i in np.arange(len(R)):
    plot_time[i] = 1.0*(i*dump_time)/1000.0 # convert fs to ps

# save R and Cv to files
RCout.write("t R\n")
if VELOCITY:
    CVout.write("t Cv\n")
for i in np.arange(len(R)):
    RCout.write("{0} {1}\n".format(plot_time[i],R[i]))
    if VELOCITY:
        CVout.write("{0} {1}\n".format(plot_time[i],Cv[i]))
    if MSD:
        MSDout.write("{0} {1}\n".format(plot_time[i],msdisp[i]))
RCout.close()
if VELOCITY:
    CVout.close()
if MSD:
    MSDout.close()
logout.close()

print "Wrote occupation time to {0}".format(RCoutputfile)
if VELOCITY:
    print "Wrote velocity correlation to {0}".format(CVoutputfile)
if MSD:
    print "Wrote mean square disp to {0}".format(MSDoutputfile)
print "Wrote log to {0}".format(logoutputfile)

#plot_rc("R_C_d{0}.png".format(z_dist))

#if VELOCITY:
    #plot_cv("C_V_d{0}.png".format(z_dist))
