//
//  correlation.c
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
#include "correlation.h"
#include <vector>
#include <string>
#include <string.h>
#include <sstream>

#define OID 2

// puts results in a pre-constructed vector

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss,item,delim)) {
        elems.push_back((item));
        //elems.push_back(std::stof(item));
    }
    //return elems;
}

bool is_in_box(double x, double y, double z, double x0, double x1, double y0, double y1, double z0, double z1) {
    return ((x > x0) && (x < x1) && (y > y0) && (y < y1) && (z > z0) && (z < z1));
}

int get_idx(int step,int particle_idx,int k,int num_particles) {
    return step*num_particles*3+particle_idx*3+k;
}

int get_idx(int step,int particle_idx,int num_particles) {
    return step*num_particles+particle_idx;
}

double sq_displacement(double x, double y, double z, double x0, double y0, double z0) {
    double sq_disp=(x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
    return sq_disp;
}

void read_traj(std::string file_name,int* theta, double* velocities, double* positions,
        double* R, double* Cv, double* MSD, double* plot_time,
        int skip, int num_particles, int num_timesteps, int max_meas_time,
        double x0, double x1, double y0, double y1, double z0, double z1, std::string suffix,
        int O_id) {
    printf("Reading %s\ntotal time %d, meas time %d, %d particles\n",
            &file_name[0],num_timesteps,max_meas_time,num_particles);
    printf("x0 x1 y0 y1 z0 z1\n%f %f %f %f %f %f\n",x0,x1,y0,y1,z0,z1);
    using namespace std;
    filebuf file_buffer;
    string str;
    if (file_buffer.open(file_name,std::ios::in))
    {

        istream file_stream(&file_buffer);
        vector<string> line_split;

        const int num_atoms=5184;

        int N_occu_total=0;
        double omega_total=0;
        int particle_idx;

        // m x n, array[i][j] = array[i*m + j]
        // l x m x n, array[i][j][k] = array[i*l*m + j*m + k]
        theta=(int*) malloc(num_particles*num_timesteps*sizeof(int));
        velocities=(double*) malloc(3*num_particles*num_timesteps*sizeof(double));
        positions=(double*) malloc(3*num_timesteps*num_particles*sizeof(double));
        R=(double*) malloc(max_meas_time*sizeof(double));
        Cv=(double*) malloc(max_meas_time*sizeof(double));
        MSD=(double*) malloc(max_meas_time*sizeof(double));
        plot_time=(double*) malloc(max_meas_time*sizeof(double));
        int current_time;

        //file_stream.seekg(0, file_stream.beg);
        for (int step=0; step<num_timesteps;++step) {
            file_stream.ignore(1000,'\n');
            for (int atom=0; atom<num_particles; ++atom) {
                getline(file_stream,str);
                //cout << str << endl;
                split(str,' ',line_split);
                // first atom is always oxygen
                {

                    particle_idx=stof(line_split[1])-1;
                    int my_2d_index=get_idx(step,particle_idx,num_particles);
                    int my_3d_index=get_idx(step,particle_idx,0,num_particles);
                    //cout << "my_3d_index: " << my_3d_index << endl;

                    for (int k=0; k<3; ++k) {
                        //int my_index=get_idx(step,particle_idx,k,num_particles);
                        //cout << "my_index: " << my_index << endl;
                        positions[my_3d_index+k] = stof(line_split[2+k]);
                        velocities[my_3d_index+k] = stof(line_split[5+k]);
                        //cout << positions[my_index] << " ";
                        //cout << positions[my_3d_index+k] << endl;

                    }
                    // check if water is in confined space
                    if (is_in_box(positions[my_3d_index],positions[my_3d_index+1],positions[my_3d_index+2],
                            x0,x1,y0,y1,z0,z1) )
                    {
                        // printf("I'm in the box, coordinates %f %f %f, my id: %d\n",
                        //         positions[my_3d_index],positions[my_3d_index+1],positions[my_3d_index+2],
                        //         particle_idx);
                        theta[my_2d_index]=1;
                        N_occu_total+=1;
                        for (int k=0; k<3; ++k) {
                            omega_total+=velocities[my_3d_index+k]*velocities[my_3d_index+k];
                        }

                    }

                }
                //cout << line_split[0] << endl;
                //cout << "---------------" << endl;
                line_split.clear();

                // last two lines are hydrogens
                file_stream.ignore(1000,'\n');
                file_stream.ignore(1000,'\n');
            }

            //printf("-------------------\nmeasuring correlation at step %d\n",step);
            // computing correlation functions
            for (int step_meas=0; step_meas < min(step,max_meas_time); ++step_meas) {
                double occu_tmp=0;
                double vel_tmp=0;
                double sq_disp_tmp=0;

                int t0=step-step_meas;
                for (int particle=0; particle < num_particles; ++particle) {
                    int my_2d_index_t0=get_idx(t0,particle,num_particles);
                    int my_2d_index=get_idx(step,particle,num_particles);
                    int my_3d_index_t0=get_idx(t0,particle,0,num_particles);
                    int my_3d_index=get_idx(step,particle,0,num_particles);
                    if (theta[my_2d_index_t0]==1) {
                        occu_tmp+=theta[my_2d_index];

                        for (int k=0; k<3; ++k) {
                            vel_tmp += velocities[my_3d_index_t0+k]*velocities[my_3d_index+k];
                        }
                        sq_disp_tmp += sq_displacement(
                                positions[my_3d_index],positions[my_3d_index+1],positions[my_3d_index+2],
                                positions[my_3d_index_t0],positions[my_3d_index_t0+1],positions[my_3d_index_t0+2]);
                    }
                }
                //printf("Found %f correlations from %d steps ago\n",occu_tmp,step_meas);
                R[step_meas]+=occu_tmp;
                Cv[step_meas]+=vel_tmp;
                MSD[step_meas]+=sq_disp_tmp;
            }

        }

        // normalizing R, Cv, and MSD
        double N_occu_ave = N_occu_total*1.0/num_timesteps;
        int t0_max;
        double omega;
        int dump_time=20;
        for (int i=0; i<max_meas_time; ++i) {
            t0_max=num_timesteps-i;
            R[i] = R[i]*1.0/N_occu_ave/(t0_max);
            omega = 1.0*omega_total/N_occu_ave/(num_timesteps);
            Cv[i] = Cv[i]*1.0/N_occu_ave/(t0_max)/omega;
            MSD[i] = MSD[i]/(t0_max);
            plot_time[i] = 1.0*(i*dump_time)/1000.0;
        }

        ofstream RCout;
        RCout.open("results/R_C_out_s"+to_string(skip)+"_m"+to_string(max_meas_time)+"_t"+to_string(num_timesteps)+"-"+suffix);
        ofstream CVout;
        CVout.open("results/C_V_out_s"+to_string(skip)+"_m"+to_string(max_meas_time)+"_t"+to_string(num_timesteps)+"-"+suffix);
        ofstream MSDout;
        MSDout.open("results/MSD_out_s"+to_string(skip)+"_m"+to_string(max_meas_time)+"_t"+to_string(num_timesteps)+"-"+suffix);
        RCout << "t R" << endl;
        CVout << "t Cv" << endl;
        MSDout << "t MSD" << endl;
        for (int i=0; i< max_meas_time; ++i) {
            RCout << plot_time[i] << " " << R[i] << endl;
            CVout << plot_time[i] << " " << Cv[i] << endl;
            MSDout << plot_time[i] << " " << MSD[i] << endl;
        }
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
    double* MSD;
    double* plot_time;
    double x0=-6.05+15;
    double x1=6.05+15;
    double y0=-5.25+15;
    double y1=5.25+15;
    double z0=16.1;
    double z1=22.90;
    int skip_steps=0;
    int num_particles=1728;
    int num_timesteps=200;
    int max_meas_time=100;
    int O_id=2;
    clock_t start, time;
    string suffix;

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
        if (strcmp("-e", argv[i])==0) {
            suffix = string(argv[i+1]);
        }
        if (strcmp("-o", argv[i])==0) {
            O_id = atoi(argv[i+1]);
        }
    }
    start=clock();

    read_traj(file_name,theta,velocities,positions,R,Cv,MSD,
            plot_time,skip_steps,num_particles,num_timesteps,max_meas_time,
            x0,x1,y0,y1,z0,z1,suffix,O_id);
    time=clock()-start;
    time=time/CLOCKS_PER_SEC;
    cout << "Took " << time << " sec to run" << endl;
    cout << "--------------------" << endl;
    return 0;
}
