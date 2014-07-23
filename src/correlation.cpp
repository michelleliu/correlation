//
//  correlation.cpp
//
//  code to calculate correlation function
//  from LAMMPS trajectory
//
//  Created by Michelle Liu (2014)
//  Chemical and Biomolecular Engineering, UC Berkeley
//

#include <fstream>
#include <iostream>
#include <cmath>
#include <array>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <string.h>
#include <sstream>
#include "common.h"

#define OID 1

// puts results in a pre-constructed vector

void read_traj(std::string file_name,int* theta, double* velocities, double* positions,
        double* R, double* Cv, double* Cv_int, double* MSD, double* Cm,
        double* mu, double* Mnet, double* plot_time, double* full_time,
        int skip, int dump_time, int num_particles, int num_timesteps, int max_meas_time,
        double x0, double x1, double y0, double y1, double z0, double z1,
        std::string suffix, int O_id, bool VERBOSE) {
    // requires that lammps trajectory be of the following format(ish):
    // dump 1 all custom 1 ${jobname}.t${jobid}.traj id mol type x y z vx vy vz
    printf("Reading %s\ntotal time %d, meas time %d, %d particles\nsuffix: %s\n",
            &file_name[0],num_timesteps,max_meas_time,num_particles,&suffix[0]);
    printf("x0 x1 y0 y1 z0 z1\n%f %f %f %f %f %f\n",x0,x1,y0,y1,z0,z1);
    using namespace std;
    filebuf file_buffer;
    string str;
    if (file_buffer.open(file_name,std::ios::in))
    {

        istream file_stream(&file_buffer);
        vector<double> line_split;

        // parse beginning lines
        file_stream.ignore(1000,'\n');
        file_stream.ignore(1000,'\n');
        file_stream.ignore(1000,'\n');
        getline(file_stream,str);
        int num_atoms=stoi(str);
        //num_atoms=5;
        int N_occu_total=0;
        double omega_total=0;
        int particle_idx;

        // m x n, array[i][j] = array[i*m + j]
        // l x m x n, array[i][j][k] = array[i*l*m + j*m + k]
        theta=(int*) malloc(num_particles*num_timesteps*sizeof(int));

        if (VERBOSE) {
            printf("Theta initial: [");
            for (int i=0; i<max_meas_time; ++i) {
                int my_2d_index=get_idx(i,1394,num_particles);
                printf("%d  %d\n",my_2d_index,theta[my_2d_index]);
            }
        }

        velocities=(double*) malloc(3*num_particles*num_timesteps*sizeof(double));
        positions=(double*) malloc(3*num_timesteps*num_particles*sizeof(double));
        R=(double*) malloc(max_meas_time*sizeof(double));
        Cv=(double*) malloc(max_meas_time*sizeof(double));
        Cv_int=(double*) malloc(max_meas_time*sizeof(double));
        MSD=(double*) malloc(max_meas_time*sizeof(double));
        plot_time=(double*) malloc(max_meas_time*sizeof(double));
        full_time=(double*) malloc(num_timesteps*sizeof(double));
        int current_time;

        // for dipole
        //mu=(double*) malloc(3*num_particles*num_timesteps*sizeof(double)); // molecular dipoles
        //Mnet=(double*) malloc(3*num_timesteps*sizeof(double)); // net dipole of confined region
        //Cm=(double*) malloc(max_meas_time*sizeof(double)); // dipole-dipole correlation
        //double q;

        // seek to beginning of file
        file_stream.seekg(0, file_stream.beg);
        for (int step=0; step<num_timesteps;++step) {
            file_stream.ignore(1000,'\n');
            getline(file_stream,str);
            current_time=stoi(str);
            file_stream.ignore(1000,'\n');
            file_stream.ignore(1000,'\n');
            file_stream.ignore(1000,'\n');
            file_stream.ignore(1000,'\n');
            file_stream.ignore(1000,'\n');
            file_stream.ignore(1000,'\n');
            file_stream.ignore(1000,'\n');

            for (int atom=0; atom<num_atoms; ++atom) {
                getline(file_stream,str);
                split(str,' ',line_split);

                // if atom is hydrogen
                // read in position and charge for dipole
                if (line_split[2]==O_id+1)
                {
                    particle_idx=line_split[1]-1; // starts at 0
                    int my_3d_index=get_idx(step,particle_idx,0,num_particles);
                    //q=line_split[9];
                    //for (int k=0; k<3; ++k) {
                    //    mu[my_3d_index+k] += q*line_split[3+k]; // xyz components of dipole
                    //}
                }

                // check atom is oxygen
                // for correlations
                else if (line_split[2]==O_id) {
                    particle_idx=line_split[1]-1;
                    int my_2d_index=get_idx(step,particle_idx,num_particles);
                    int my_3d_index=get_idx(step,particle_idx,0,num_particles);

                    //q=line_split[9]; // dipole

                    for (int k=0; k<3; ++k) {
                        positions[my_3d_index+k] = line_split[3+k];
                        velocities[my_3d_index+k] = line_split[6+k];
                        //mu[my_3d_index+k] += q*line_split[3+k]; // xyz components of dipole

                    }
                    // check if water is in confined space
                    if (check_in_box(positions[my_3d_index],positions[my_3d_index+1],positions[my_3d_index+2],
                            x0,x1,y0,y1,z0,z1) )
                    {
                        if (VERBOSE) {
                            printf("I'm in the box at step %d, "
                                    "coordinates %f %f %f, "
                                    "my id: %d, my_2d_idx: %d\n",
                                    step,positions[my_3d_index],positions[my_3d_index+1],
                                    positions[my_3d_index+2],particle_idx,my_2d_index);
                        }
                        theta[my_2d_index]=1;
                        N_occu_total+=1;
                        for (int k=0; k<3; ++k) {
                            omega_total+=velocities[my_3d_index+k]*velocities[my_3d_index+k];
                        }

                    }

                }
                line_split.clear();
            }

            // calculating net dipole for confined region
            //for (int particle=0; particle < num_particles; ++particle) {
            //        int my_2d_index=get_idx(step,particle,num_particles);
            //        int my_3d_index=get_idx(step,particle,0,num_particles);
            //        if (theta[my_2d_index]==1) {
            //            for (int k=0; k<3; ++k) {
            //                Mnet[3*step+k] += mu[my_3d_index+k];
            //            }
            //        }

            //}

            if (VERBOSE) {
                printf("-------------------\nmeasuring correlation at step %d\n",step);
            }
            // computing correlation functions
            for (int step_meas=0; step_meas < min(step,max_meas_time); ++step_meas) {
                int occu_tmp=0;
                double vel_tmp=0;
                double sq_disp_tmp=0;

                int t0=step-step_meas;
                //for (int k=0; k<3; ++k) {
                //    Cm[step_meas]+=Mnet[3*step+k]*Mnet[3*t0+k];
                //}
                for (int particle=0; particle < num_particles; ++particle) {
                    // initial and current 2-D times
                    int my_2d_index_t0=get_idx(t0,particle,num_particles);
                    int my_2d_index=get_idx(step,particle,num_particles);
                    // initial and current 3-D times
                    int my_3d_index_t0=get_idx(t0,particle,0,num_particles);
                    int my_3d_index=get_idx(step,particle,0,num_particles);
                    // check if particle in box at initial time
                    if (theta[my_2d_index_t0]==1) {
                        // add 1 to occu_tmp if particle in box at current time
                        occu_tmp+=theta[my_2d_index];

                        for (int k=0; k<3; ++k) {
                            vel_tmp += velocities[my_3d_index_t0+k]*velocities[my_3d_index+k];
                        }
                        sq_disp_tmp += sq_displacement(
                                positions[my_3d_index],positions[my_3d_index+1],positions[my_3d_index+2],
                                positions[my_3d_index_t0],positions[my_3d_index_t0+1],positions[my_3d_index_t0+2]);
                    }
                }
                if (VERBOSE) {
                    printf("Found %f correlations from %d steps ago\n",occu_tmp,step_meas);
                }
                R[step_meas]+=occu_tmp;
                Cv[step_meas]+=vel_tmp;
                MSD[step_meas]+=sq_disp_tmp;
            }

        }

        // normalizing R, Cv, and MSD
        int t0_max;
        double omega;
        printf("Pre-normalization R: [");
        for (int i=0; i<29; ++i) {
            printf("%f, ",R[i]);
        }
        printf("%f]\n",R[30]);

        if (VERBOSE) {
            printf("Theta final: [");
            for (int i=0; i<max_meas_time; ++i) {
                int my_2d_index=get_idx(i,1394,num_particles);
                printf("%d  %d\n",my_2d_index,theta[my_2d_index]);
            }
        }

        printf("N_occu_total = %d\n",N_occu_total,num_timesteps);
        double N_occu_ave = 1.0*N_occu_total/num_timesteps;
        omega = (1.0*omega_total/N_occu_ave)/(num_timesteps);
        printf("N_occu_ave = %f\n",N_occu_ave);
        double R0=R[0];
        for (int i=0; i<max_meas_time; ++i) {
            t0_max=num_timesteps-i;
            if (VERBOSE) {
                printf("%d\n",t0_max);
            }
            R[i] = (1.0*R[i]/N_occu_ave)/(t0_max);
            //R[i] = R[i]/R0; // alternative normalization
            //can't use because of truncated time window
            Cv_int[i] = Cv[i]/t0_max/N_occu_ave;
            Cv[i] = ((Cv[i]/N_occu_ave)/(t0_max))/omega;
            MSD[i] = (MSD[i]/N_occu_ave)/(t0_max);
            plot_time[i] = (i*dump_time)/1000.0; // now in ps
        }
        for (int i=0; i<num_timesteps; ++i) {
            full_time[i] = (i*dump_time)/1000.0; // now in ps
        }

        double area;
        integrate(Cv_int, plot_time, 0, max_meas_time, &area);
        double D;
        D=area*100/3;
        printf("D = %f cm^2/s\n",D);


        ofstream RCout;
        RCout.open("results/R_C-"+suffix+"-w"+to_string(num_particles)+"_m"+to_string(max_meas_time)
                +"_t"+to_string(num_timesteps)+"_d"+to_string(dump_time)+"fs");
        ofstream CVout;
        CVout.open("results/C_V-"+suffix+"-w"+to_string(num_particles)+"_m"+to_string(max_meas_time)
                +"_t"+to_string(num_timesteps)+"_d"+to_string(dump_time)+"fs");
        ofstream MSDout;
        MSDout.open("results/MSD-"+suffix+"-w"+to_string(num_particles)+"_m"+to_string(max_meas_time)
                +"_t"+to_string(num_timesteps)+"_d"+to_string(dump_time)+"fs");
        //ofstream CMout;
        //CMout.open("results/C_M-"+suffix+"-w"+to_string(num_particles)+"_m"+to_string(max_meas_time)
        //        +"_t"+to_string(num_timesteps)+"_d"+to_string(dump_time)+"fs");
        //ofstream MNETout;
        //MNETout.open("results/Mnet-"+suffix+"-w"+to_string(num_particles)+"_m"+to_string(max_meas_time)
        //        +"_t"+to_string(num_timesteps)+"_d"+to_string(dump_time)+"fs");
        RCout << "#t R" << endl;
        CVout << "#t Cv" << endl;
        MSDout << "#t MSD" << endl;
        //CMout << "#t Cm" << endl;
        for (int i=0; i< max_meas_time; ++i) {
            RCout << plot_time[i] << " " << R[i] << endl;
            CVout << plot_time[i] << " " << Cv[i] << endl;
            MSDout << plot_time[i] << " " << MSD[i] << endl;
            //CMout << plot_time[i] << " " << Cm[i] << endl;
        }
        //for (int i=0; i< num_timesteps; ++i) {
        //    MNETout << full_time[i] << " " << Mnet[3*i+0] << " " << Mnet[3*i+1] << " " << Mnet[3*i+2] << endl;
        //}
        RCout.close();
        CVout.close();
        MSDout.close();
        //CMout.close();
        //MNETout.close();
        file_buffer.close();
    }
    else {
        cout << "ERROR opening file!" << endl;
    }
}

int main(int argc, const char * argv[]) {
    using namespace std;
    int* theta;
    double* velocities;
    double* positions;
    double* R;
    double* Cv;
    double* Cv_int;
    double* MSD;
    double* Cm;
    double* Mnet;
    double* mu;
    double* plot_time;
    double* full_time;
    double x0=-6.05;
    double x1=6.05;
    double y0=-5.25;
    double y1=5.25;
    double z0=16.1;
    double z1=22.90;
    int skip_steps=0;
    int num_particles=1000;
    int num_timesteps=200;
    int max_meas_time=100;
    int dump_time=20;
    int O_id=2;
    clock_t start, time;
    string suffix;
    bool VERBOSE;

    string file_name=argv[argc-1];
    for (int i = 0; i < (argc-1); ++i) {
        if (strcmp("-t", argv[i])==0) {
            num_timesteps = atoi(argv[i+1]);
        }
        if (strcmp("-d", argv[i])==0) {
            double z_dist = atof(argv[i+1]);
            double z_center = 19.5;
            z0=z_center-z_dist/2;
            z1=z_center+z_dist/2;
        }
        if (strcmp("-m", argv[i])==0) {
            max_meas_time = atoi(argv[i+1]);
        }
        if (strcmp("-s", argv[i])==0) {
            skip_steps = atoi(argv[i+1]);
        }
        if (strcmp("-r", argv[i])==0) {
            dump_time = atoi(argv[i+1]);
        }
        if (strcmp("-e", argv[i])==0) {
            suffix = string(argv[i+1]);
        }
        if (strcmp("-o", argv[i])==0) {
            O_id = atoi(argv[i+1]);
        }
        if (strcmp("-w", argv[i])==0) {
            num_particles = atoi(argv[i+1]);
        }
        if (strcmp("-C", argv[i])==0) {
            x0+=15;
            x1+=15;
            y0+=15;
            y1+=15;
        }
        if (strcmp("-V", argv[i])==0) {
            VERBOSE=1;
        }
        else { VERBOSE=0; }
    }
    start=clock();

    read_traj(file_name,theta,velocities,positions,
            R,Cv,Cv_int,MSD,Cm,mu,Mnet,
            plot_time,full_time,skip_steps,dump_time,
            num_particles,num_timesteps,max_meas_time,
            x0,x1,y0,y1,z0,z1,suffix,O_id,VERBOSE);
    time=clock()-start;
    time=time/CLOCKS_PER_SEC;
    cout << "Took " << time << " sec to run" << endl;
    cout << "--------------------" << endl;
    return 0;
}
