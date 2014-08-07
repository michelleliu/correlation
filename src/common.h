#ifndef __CORR_COMMON_H__
#define __CORR_COMMON_H__

#include <vector>

void split(const std::string &s, char delim, std::vector<double> &elems) ;
bool is_in_box(double x, double y, double z, double x0, double x1,
        double y0, double y1, double z0, double z1) ;
bool check_in_box(double x, double y, double z, double x0, double x1,
        double y0, double y1, double z0, double z1) ;
double wrap_pbc(double xu, double xlo, double xhi) ;
int get_idx(int step,int particle_idx,int atom_idx,int k,int num_particles) ;
int get_idx(int step,int particle_idx,int k,int num_particles) ;
int get_idx(int step,int particle_idx,int num_particles) ;
double distance(double x, double y, double z, double x0, double y0, double z0) ;
double distance_pbc(double x, double y, double z, double x0, double y0, double z0, double x_width, double y_width, double z_width) ;
double sq_displacement(double x, double y, double z, double x0, double y0,
        double z0) ;
double angle(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2) ;
void integrate(double* func, double* x, int start, int end, double* area) ;

#endif
