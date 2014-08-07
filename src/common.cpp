#include <fstream>
#include <iostream>
#include <cmath>
#include <array>
#include <time.h>
#include <assert.h>
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <string>
#include <string.h>
#include <sstream>
#include "common.h"
#include <math.h>
#define PI 3.141592654

// puts results in a pre-constructed vector
void split(const std::string &s, char delim, std::vector<double> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss,item,delim)) {
        elems.push_back(std::stof(item));
    }
    //return elems;
}

bool is_in_box(double x, double y, double z, double x0, double x1, double y0, double y1, double z0, double z1) {
    return 1;
}

bool check_in_box(double x, double y, double z, double x0, double x1, double y0, double y1, double z0, double z1) {
    return ((x > x0) && (x < x1) && (y > y0) && (y < y1) && (z > z0) && (z < z1));
}

double wrap_pbc(double xu, double xlo, double xhi) {
    double xwidth = xhi - xlo;
    double x = xlo + xu - xwidth * floor( (xu - xlo) / xwidth );
    return x;
}
// hydrogens index (for h_positions)
int get_idx(int step,int particle_idx,int atom_idx,int k,int num_particles) {
    int j=atom_idx%3-1;
    //printf("atom_idx=%d, j=%d\n",atom_idx,j);
    assert(j==0 || j==1);
    return 2*step*num_particles*3+(2*particle_idx+j)*3+k;
}

// 3d index (for velocities, positions)
int get_idx(int step,int particle_idx,int k,int num_particles) {
    return step*num_particles*3+particle_idx*3+k;
}

// 2d index
int get_idx(int step,int particle_idx,int num_particles) {
    return step*num_particles+particle_idx;
}

double sq_displacement(double x, double y, double z, double x0, double y0, double z0) {
    double sq_disp=(x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
    return sq_disp;
}

double distance(double x, double y, double z, double x0, double y0, double z0) {
    double sq_disp=(x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
    return sqrt(sq_disp);
}

double distance_pbc(double x, double y, double z, double x0, double y0, double z0, double x_width, double y_width, double z_width) {
    //printf("%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",x,y,z,z0,y0,z0,x_width,y_width,z_width);
    double del_x = (x-x0);
    double del_y = (y-y0);
    double del_z = (z-z0);
    if (del_x < -x_width*0.5) {del_x = del_x + x_width;}
    if (del_x >= x_width*0.5) {del_x = del_x - x_width;}
    if (del_y < -y_width*0.5) {del_y = del_y + y_width;}
    if (del_y >= y_width*0.5) {del_y = del_y - y_width;}
    if (del_z < -z_width*0.5) {del_z = del_z + z_width;}
    if (del_z >= z_width*0.5) {del_z = del_z - z_width;}
    double sq_disp=del_x*del_x+del_y*del_y+del_z*del_z;
    return sqrt(sq_disp);
}

double angle(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2) {
    double a [3];
    double b [3];
    a[0]=x0-x1;
    a[1]=y0-y1;
    a[2]=z0-z1;
    b[0]=x2-x1;
    b[1]=y2-y1;
    b[2]=z2-z1;
    double dotprod = 0;
    for (int i=0; i<3; ++i) {
        dotprod += a[i]*b[i];
    }
    double a_mag = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    double b_mag = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
    double phi = acos( dotprod / (a_mag*b_mag) ) * 180.0 / PI;
    return phi;
}

void integrate(double* func, double* x, int start, int end, double* area) {
    *area = 0;
    double height = 0;
    height += (func[start] + func[end-1])/2;
    for (int i = start+1; i < end-1; i ++) {
        height += func[i];
    }
    *area = height*(x[start+1]-x[start]);
}

