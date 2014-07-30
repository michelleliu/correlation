//
//  analysis.cpp
//
//  main code for correlation functions and trajectory analysis
//  correlation, dipole, hbonds, dielectric, compressibility
//  from LAMMPS trajectory of water (confined and bulk)
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

int main(int argc, const char * argv[]) {
    using namespace std;
    string file_name=argv[argc-1];

    int* theta;
    double* velocities;
    double* positions;

    double* R;
    double* Cv;
    double* Cv_int;
    double* MSD;

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
