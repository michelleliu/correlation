//
//  test.cpp
//
//  testing code
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

int get_idx(int step,int particle_idx,int k,int num_particles, int max_meas_time) {
    step=step%max_meas_time;
    return step*num_particles*3+particle_idx*3+k;
}

int get_idx(int step,int particle_idx,int num_particles, int max_meas_time) {
    step=step%max_meas_time;
    return step*num_particles+particle_idx;
}

int my_3d_idx, my_2d_idx;
int main(int argc, const char * argv[]) {
    int max_meas_time=4;
    int num_particles=3;
    int num_timesteps=10;
    printf("%d steps, %d meas_time, %d particles\n",num_timesteps,max_meas_time,num_particles);
    printf("3d matrix %d x %d x %d\n",max_meas_time,num_particles,3);
    printf("2d matrix %d x %d\n...................\n",max_meas_time,num_particles);
    for (int step=0; step<num_timesteps; ++step) {
        for (int particle_idx=0;particle_idx<num_particles;++particle_idx) {
            for (int k=0;k<3;++k) {
                int my_3d_idx=get_idx(step, particle_idx,k,num_particles,max_meas_time);
                printf("step %d, particle_idx %d, x[%d], my_3d_idx: %d\n",step,particle_idx,k,my_3d_idx);
            }
            int my_2d_idx=get_idx(step, particle_idx,num_particles,max_meas_time);
            printf("step %d, particle_idx %d, my_2d_idx: %d\n----------\n",step,particle_idx,my_2d_idx);
        }
        printf("=====================================\n");
    }
}
