//
//  hbonds.cpp
//
//  code to analyze hbonds
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
#include <stdio.h>
#include <vector>
#include <string>
#include <string.h>
#include <sstream>
#include <assert.h>
#include <tuple>
#include "common.h"
#include <algorithm>

#include <functional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#define OID 1
#define BIN_SIZE 3*3.6

struct pairhash{
    size_t operator()(const std::pair<int, int> &p) const {
        return
            std::hash<unsigned int>()(p.first) ^
            std::hash<unsigned int>()(p.second);
    }
};
typedef std::unordered_map<std::pair<int, int>, std::vector<int>*, pairhash> hbond_map;
//typedef std::unordered_map<int, std::vector<int>*> hbond_int_map;
typedef std::unordered_map<std::pair<int, int>, std::vector<int>*, pairhash>::iterator hbond_map_itr;
//typedef std::unordered_map<int, std::vector<int>*>::iterator hbond_int_map_itr;

// puts results in a pre-constructed vector
int get_bin_idx(std::tuple<int, int, int> block_id, const int mb, const int nb) {
    using namespace std;
    return (get<0>(block_id)) * (mb) * (nb) + (get<1>(block_id)) * (nb) + (get<2>(block_id)) ;
}
std::tuple<int, int, int> get_block_id(int bin_idx, const int mb, const int nb) {
    int z_block = ( bin_idx ) % nb;
    //assert((( bin_idx ) - z_block) % nb==0);
    int y_block = ( (( bin_idx ) - z_block) / nb ) % mb;
    //assert(( bin_idx - z_block - nb*y_block) % (mb * nb)==0);
    int x_block = ( bin_idx - z_block - nb*y_block) / (mb * nb);
    return std::make_tuple(x_block,y_block,z_block);
}
int get_O_mol_idx(int atom_idx) {
    assert ((atom_idx) % 3 == 0);
    return (atom_idx) / 3;
}
int get_O_atom_idx(int mol_idx) {
    return (mol_idx) * 3;
}
int get_H_atom_idx(int mol_idx, int i) {
    // return atom index of hydrogen given mol_idx and i (1 or 2)
    return (mol_idx) * 3 + i;
}
void bin_particle(std::vector<int>** bins, int particle_idx,
        double particle_x, double particle_y, double particle_z,
        int mb, int nb, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    //printf("begin binning\n");
    using namespace std;
    tuple<int, int, int> block_id;
    int* block_x=&get<0>(block_id);
    int* block_y=&get<1>(block_id);
    int* block_z=&get<2>(block_id);
    double bin_size_x=(xmax-xmin)/mb;
    double bin_size_y=(ymax-ymin)/mb;
    double bin_size_z=(zmax-zmin)/nb;
    //printf("bin sizes: %f %f %f\n",bin_size_x,bin_size_y,bin_size_z);
    *block_x = (int) floor((particle_x-xmin) / (bin_size_x));
    *block_y = (int) floor((particle_y-ymin) / (bin_size_y));
    *block_z = (int) floor((particle_z-zmin) / (bin_size_z));
    if (*block_x < 0) {*block_x = 0;}
    if (*block_y < 0) {*block_y = 0;}
    if (*block_z < 0) {*block_z = 0;}
    if (*block_x >=mb) {*block_x = mb-1;}
    if (*block_y >=mb) {*block_y = mb-1;}
    if (*block_z >=nb) {*block_z = nb-1;}

    int bin_idx = get_bin_idx(block_id, mb, nb);
    //printf("bin_idx: %d, block %d %d %d\n",bin_idx, get<0>(block_id),get<1>(block_id),get<2>(block_id));
    bins[bin_idx]->push_back(particle_idx);
    // printf("adding %d to bin %d, block %d %d %d, position: %f, %f, %f\n",
    //         particle_idx,bin_idx,
    //         *block_x,*block_y,*block_z,
    //         particle_x,particle_y,particle_z);
}
void clear_bins(std::vector<int>** bins, const int mb, const int nb) {
    for (int i = 0; i < (nb * mb * mb); ++i) {
        //omp_set_lock(&locks[i]);
        bins[i]->clear();
        //omp_unset_lock(&locks[i]);
    }
}
void clear_nbrs(std::vector<int>** nbrs, const int num_particles) {
    for (int i = 0; i < (num_particles); ++i) {
        //omp_set_lock(&locks[i]);
        nbrs[i]->clear();
        //omp_unset_lock(&locks[i]);
    }
}
// using "block" and "bin" to mean the same thing, only "blocks" have (x,y,z) names
void get_neighbors(const std::tuple<int, int, int>& block_id, std::vector<std::tuple<int, int, int> >* neighbors, int mb, int nb) {
    using namespace std;
    //printf("For block (%d, %d, %d) nbrs: ", get<0>(block_id), get<1>(block_id), get<2>(block_id));
    for (int i = get<0>(block_id) - 1; i <= get<0>(block_id) + 1; ++i) {
        for (int j = get<1>(block_id) - 1; j <= get<1>(block_id) + 1; ++j) {
            for (int k = get<2>(block_id) - 1; k <= get<2>(block_id) + 1; ++k) {
                if (!(i == get<0>(block_id) && j == get<1>(block_id) && k == get<2>(block_id) )) {
                    int ii=i;
                    int jj=j;
                    int kk=k;
                    if (ii < 0) {ii = ii+mb;}
                    if (jj < 0) {jj = jj+mb;}
                    if (kk < 0) {kk = kk+nb;}
                    if (ii >=mb) {ii = ii-mb;}
                    if (jj >=mb) {jj = jj-mb;}
                    if (kk >=nb) {kk = kk-nb;}
                    //if ( i >= 0 && j >= 0 && k >= 0 ) {
                    //    if ( i < mb && j < mb && k < nb ) {
                            //printf(" (%d, %d, %d), ", i, j, k);
                    neighbors->push_back(make_tuple(ii, jj, kk));
                    //    }
                    //}
                }
            }
        }
    }
    //printf("\n");
}

bool is_H_bond(int donor_O, int donor_H, int acceptor_O, int step, int num_particles,
        double oh_cutoff, double phi_cutoff, double* positions, double* h_positions,
        double x_width, double y_width, double z_width) {

    //printf("testing hbond between %d, %d, %d, ",donor_O,donor_H,acceptor_O,step);
    int particle_idx=get_O_mol_idx(donor_O);
    int H_atom_idx=donor_O+donor_H;

    int donor_3d_index = get_idx(step,particle_idx,0,num_particles);
    int donor_H_3d_index = get_idx(step,particle_idx,H_atom_idx,0,num_particles);
    int acceptor_3d_index = get_idx(step,get_O_mol_idx(acceptor_O),0,num_particles);

    double oh_distance=distance_pbc(
            positions[acceptor_3d_index],positions[acceptor_3d_index+1],positions[acceptor_3d_index+2], // acceptor O
            h_positions[donor_H_3d_index],h_positions[donor_H_3d_index+1],h_positions[donor_H_3d_index+2],
            x_width, y_width, z_width) ; // donor H
    //printf("oh_distance: %f, ",oh_distance);
    if (oh_distance < oh_cutoff) {
        double phi = angle(
                positions[donor_3d_index],positions[donor_3d_index+1],positions[donor_3d_index+2] , // donor O
                h_positions[donor_H_3d_index],h_positions[donor_H_3d_index+1],h_positions[donor_H_3d_index+2] , // donor H
                positions[acceptor_3d_index],positions[acceptor_3d_index+1],positions[acceptor_3d_index+2] ); // acceptor O
        //printf("phi: %f, ",phi);
        if (phi > phi_cutoff) {
            //printf("yes hbond\n");
            return 1;
        } else {
            //printf("no hbond\n");
            return 0;
        }
    }
    else {
        //printf("no hbond\n");
        return 0;
    }
}

void read_hbonds(std::string file_name,int* theta, double* Ch, double* positions, double* h_positions,
        std::vector<int>** bins, int mb, int nb,
        double oo_cutoff, double oh_cutoff, double phi_cutoff, hbond_map* hbonds,
        double* plot_time, double* full_time,
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
        //double omega_total=0;
        int particle_idx;
        int atom_idx;

        // m x n, array[i][j] = array[i*m + j]
        // l x m x n, array[i][j][k] = array[i*l*m + j*m + k]
        theta=(int*) malloc(num_particles*num_timesteps*sizeof(int));
        Ch=(double*) malloc(max_meas_time*sizeof(double));

        if (VERBOSE) {
            printf("Theta initial: [");
            for (int i=0; i<max_meas_time; ++i) {
                int my_2d_index=get_idx(i,1394,num_particles);
                printf("%d  %d",my_2d_index,theta[my_2d_index]);
                printf("]\n");
            }
        }

        positions=(double*) malloc(3*num_timesteps*num_particles*sizeof(double));
        h_positions=(double*) malloc(3*2*num_timesteps*num_particles*sizeof(double));
        plot_time=(double*) malloc(max_meas_time*sizeof(double));
        full_time=(double*) malloc(num_timesteps*sizeof(double));
        int current_time;
        vector<int>** nbr_particle_list;
        nbr_particle_list = new vector<int>*[(num_particles)];
        printf("%d num_particles\n",num_particles);
        for (int i=0; i<num_particles; ++i) {
            nbr_particle_list[i] = new vector<int>();
        }
        vector<int> num_hbonds(num_timesteps);
        for (int i=0; i<num_timesteps; ++i) {
            num_hbonds[i]=0;
        }
        //printf("hbonds at start:\n");
        //for (int i=0; i<num_timesteps; ++i) {
        //    printf("step %d, %d hbonds\n",i,num_hbonds[i]);
        //}

        // seek to beginning of file
        file_stream.seekg(0, file_stream.beg);
        double xmin, xmax, ymin, ymax, zmin, zmax;
        int hbonds_this_step;
        for (int step=0; step<num_timesteps;++step) {
            printf("-----------\n");
            printf("Beginning step %d\n",step);
            file_stream.ignore(1000,'\n');
            getline(file_stream,str);
            current_time=stoi(str);
            printf("TIMESTEP: %d\n",current_time);
            file_stream.ignore(1000,'\n');
            file_stream.ignore(1000,'\n');
            file_stream.ignore(1000,'\n');
            //file_stream.ignore(1000,'\n');
            //file_stream.ignore(1000,'\n');
            //file_stream.ignore(1000,'\n');
            getline(file_stream,str);
            split(str,' ',line_split);
            xmin=line_split[0];
            xmax=line_split[1];
            getline(file_stream,str);
            split(str,' ',line_split);
            ymin=line_split[0];
            ymax=line_split[1];
            getline(file_stream,str);
            split(str,' ',line_split);
            zmin=line_split[0];
            zmax=line_split[1];
            file_stream.ignore(1000,'\n');

            double x_width=xmax-xmin;
            double y_width=ymax-ymin;
            double z_width=zmax-zmin;
            printf("x,y,z dimensions: %f %f %f\n",x_width,y_width,z_width);
            hbonds_this_step=0;

            for (int atom=0; atom<num_atoms; ++atom) {
                getline(file_stream,str);
                split(str,' ',line_split);

                // if atom is hydrogen
                if (line_split[2]==O_id+1) {
                    particle_idx=line_split[1]-1; // starts at 0
                    atom_idx=line_split[0]-1; // starts at 0
                    int my_h_index=get_idx(step,particle_idx,atom_idx,0,num_particles);
                    //int my_2d_index=get_idx(step,particle_idx,num_particles);
                    //int my_3d_index=get_idx(step,particle_idx,0,num_particles);
                    for (int k=0; k<3; ++k) {
                        h_positions[my_h_index+k] = line_split[3+k];
                    }
                    //printf("Got here, hydrogen, atom %d, step %d\n", atom,step);
                }

                // check atom is oxygen
                // for correlations
                else if (line_split[2]==O_id) {
                    particle_idx=line_split[1]-1;
                    int my_2d_index=get_idx(step,particle_idx,num_particles);
                    int my_3d_index=get_idx(step,particle_idx,0,num_particles);

                    //printf("Got here, oxygen, atom: %d, my particle_idx: %d, my_3d_index: %d, step %d\n", atom, particle_idx, my_3d_index,step);
                    for (int k=0; k<3; ++k) {
                        //printf("beginning of loop\n");
                        positions[my_3d_index+k] = line_split[3+k];
                    }
                    //printf("Got here, oxygen, after storing positions\n");
                    // bin particle
                    //printf("binning particle_idx %d\n", particle_idx);
                    bin_particle(bins, particle_idx, positions[my_3d_index+0],
                            positions[my_3d_index+1], positions[my_3d_index+2],
                            mb, nb, xmin, xmax, ymin, ymax, zmin, zmax);
                    //printf("after binning particles\n");
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

                    }

                    //printf("end of oxygen\n");
                }
                line_split.clear();
            }
            printf("finished reading atoms\n");
            printf("finished binning particles\n");
            //printf("bin 5 has %d waters:\n",bins[5]->size());
            //for (int i=0; i<bins[5]->size(); ++i) {
            //    printf("%d\n", bins[5]->at(i));
            //}


            if (VERBOSE) {
                printf("-------------------\nmeasuring correlation at step %d\n",step);
            }


            // check whether iterated through this pair of bins yet
            vector<bool> check_pair(mb*mb*nb*mb*mb*nb);
            for (int i = 0; i < (mb*mb*nb*mb*mb*nb); ++i) {
                check_pair[i]=0;
                //printf("%d, ",check_pair[i]);
            }
            //printf("\n");
            printf("finished initializing check_pair\n");

            // iterate through bins to find neighbors
            printf("there are %d bins\n",nb*mb*mb);
            for (int i = 0; i < (nb*mb*mb) ; ++i) { // for all bins (outermost)
                vector<int>* this_bin = bins[i];
                vector<tuple<int, int, int> > neighbors;
                get_neighbors(get_block_id(i, mb, nb) ,  &neighbors, mb, nb); // get neighbor bins

                //tuple<int, int, int> my_block_id = get_block_id(i, mb, nb);
                //printf("got to bin %d, block %d %d %d, %d waters, %d neighbors\n",i,
                //        get<0>(my_block_id), get<1>(my_block_id), get<2>(my_block_id), this_bin->size(), neighbors.size());

                for (int k = 0; k < neighbors.size(); ++k) { // iterate over neighbor bins
                    //printf("nbr: %d\n",k);
                    int neigh_bin_idx = get_bin_idx(neighbors[k], mb, nb);
                    //printf("nbr %d, neigh_bin_idx: %d, neigh_block_id %d %d %d\n", k,
                    //        neigh_bin_idx, get<0>(neighbors[k]), get<1>(neighbors[k]), get<2>(neighbors[k]));
                    vector<int>* neighbor_bin = bins[neigh_bin_idx];

                    //printf("check_pair %d and %d? %d\n",i, neigh_bin_idx, check_pair[i*(mb*mb*nb)+neigh_bin_idx]);

                    // if this pair of bins has not been checked
                    if ((check_pair[i*(mb*mb*nb)+neigh_bin_idx])==0) {
                        // mark pair of bins as checked
                        check_pair[i*(mb*mb*nb)+neigh_bin_idx]=1;
                        check_pair[neigh_bin_idx*(mb*mb*nb)+i]=1;

                        for (int jj = 0; jj < this_bin->size(); ++jj) { // iterate over particles in this bin

                            //printf("water: %d\n",jj);
                            int my_3d_index=get_idx(step,this_bin->at(jj),0,num_particles);

                            for (int ii = 0; ii < neighbor_bin->size(); ++ii) { // iterate over particles in neighbor bin
                                //printf("nbr wat: %d\n",ii);
                                int neighbor_3d_index=get_idx(step,neighbor_bin->at(ii),0,num_particles);

                                // distance between my_particle and neigh_part
                                double r_oioj = distance_pbc(
                                        positions[my_3d_index],positions[my_3d_index+1],positions[my_3d_index+2],
                                        positions[neighbor_3d_index],positions[neighbor_3d_index+1],positions[neighbor_3d_index+2],
                                        x_width, y_width, z_width );
                                if (r_oioj < oo_cutoff) {
                                    //printf("%d and %d are neighbors\n",this_bin->at(jj),neighbor_bin->at(ii));
                                    nbr_particle_list[this_bin->at(jj)]->push_back(neighbor_bin->at(ii)); // FIXME to add simultaneously to both nbr lists (fix nbr bins)
                                    nbr_particle_list[neighbor_bin->at(ii)]->push_back(this_bin->at(jj));
                                }
                            }
                        }
                    }
                }
                neighbors.clear();
            }
            printf("finished building neighbor lists\n");

            printf("particle 10 has %d neighbors\n",nbr_particle_list[10]->size());
            printf("particle %d has %d neighbors\n",num_particles/2,nbr_particle_list[num_particles/2]->size());
            printf("particle %d has %d neighbors\n",num_particles-1,nbr_particle_list[num_particles-1]->size());
            //printf("particle 863 has %d neighbors\n",nbr_particle_list[863]->size());
            //printf("particle 1727 has %d neighbors\n",nbr_particle_list[1727]->size());
            //for (int i=0; i<nbr_particle_list[0]->size(); ++i) {
            //    printf("%d\n", nbr_particle_list[0]->at(i));
            //}

            // now have complete neighbor lists

            // analyzing timestep for H bonds
            for (int i=0; i<num_particles; ++i) {
                // find neighbors within oo_cutoff
                int donor_O = get_O_atom_idx(i);
                vector<int>* my_nbr_particles = nbr_particle_list[i];
                for (int j=0; j<my_nbr_particles->size(); ++j) {
                    int acceptor_O = get_O_atom_idx(my_nbr_particles->at(j));
                    for (int k=1; k<3; ++k) { // for each H
                        // find O in neighbors within oh_cutoff and phi_cutoff
                        if (is_H_bond(donor_O, k, acceptor_O, step, num_particles, oh_cutoff, phi_cutoff,positions,h_positions,x_width,y_width,z_width)) {
                            //printf("found a hydrogen bond between %d and %d!\n",donor_O+k,acceptor_O);
                            num_hbonds[step]+=1;
                            hbonds_this_step+=1;
                            int donor_H = donor_O+k;
                            pair<int, int> hbond_id(donor_H,acceptor_O);
                            //printf("hbonds %p, hbond_id %p, hbond_id %d %d\n",hbonds, &hbond_id, hbond_id.first, hbond_id.second);
                            hbond_map_itr itr = hbonds->find(hbond_id);
                            //printf("after finding itr\n");
                            if (itr != hbonds->end()) {
                                (itr->second)->push_back(step);
                                //printf("existing hbond at step %d\n",step);
                            } else {
                                vector<int>* nvec = new vector<int>();
                                //printf("new hbond at step %d\n",step);
                                hbonds->insert(make_pair(hbond_id,nvec));
                                //printf("got here\n");
                                //printf("size of new vector: %d\n",nvec->size());
                                nvec->push_back(step);
                                //printf("size of new vector: %d\n",nvec->size());
                            }
                            //printf("pushed back\n");
                        }
                    }
                }
            }

            // computing correlation functions
            for (int step_meas=0; step_meas < min(step,max_meas_time); ++step_meas) {
                int occu_tmp=0;

                int t0=step-step_meas;
                for (int particle=0; particle < num_particles; ++particle) {
                    // initial and current 2-D times
                    int my_2d_index_t0=get_idx(t0,particle,num_particles);
                    int my_2d_index=get_idx(step,particle,num_particles);
                    // initial and current 3-D times
                    //int my_3d_index_t0=get_idx(t0,particle,0,num_particles);
                    //int my_3d_index=get_idx(step,particle,0,num_particles);
                    // check if particle in box at initial time
                    if (theta[my_2d_index_t0]==1) {
                        // add 1 to occu_tmp if particle in box at current time
                        occu_tmp+=theta[my_2d_index];

                    }
                }
                if (VERBOSE) {
                    printf("Found %d correlations from %d steps ago\n",occu_tmp,step_meas);
                }
            }

            printf("Got to end of step %d, total %d hbonds\n",step,hbonds_this_step);
            for (int i=0; i<num_particles; ++i) {
                nbr_particle_list[i]->clear();
            }
            clear_bins(bins, mb, nb);
        }

        // calculate hbond lifetimes and correlations
        vector<int>* hb_lifetimes = new vector<int>();
        hbond_map_itr itr = hbonds->begin();
        for (int i=0;i<(hbonds)->size(); ++i) {

            vector<int>* this_hbond = (itr->second);
            int previous=-20;
            int current=-10;
            int lifetime=0;

            //for (int j=0; j<this_hbond->size(); ++j) {
                //printf("%d, ",this_hbond->at(j));
            //}
            //printf("\n");
            for (int step_meas=0; step_meas<max_meas_time; ++step_meas) {
                //printf("correlation for step_meas %d\n",step_meas);
                for (int j=0; j<this_hbond->size(); ++j) {
                    current=this_hbond->at(j);
                    assert(current>=0);
                    if ( (current-1) == previous ) {
                        lifetime+=1;
                    } else {
                        hb_lifetimes->push_back(lifetime);
                        lifetime=0;
                    }
                    int t0=current-step_meas;
                    //vector<int>::iterator it;
                    //it = find(this_hbond->begin(),this_hbond->end(),t0);
                    if (binary_search(this_hbond->begin(),this_hbond->end(),t0)) {
                        Ch[step_meas]+=1;
                        //printf("found correlation\n");
                    }
                    previous=current;
                    current=-10;
                }
            }
            itr = next(itr);
        }

        // diagnostics

        int max_lifetime=0;
        int min_lifetime=100000;

        double sum=0;
        for (int i=0; i<hb_lifetimes->size(); ++i) {
            if (hb_lifetimes->at(i)>max_lifetime) {
                max_lifetime = hb_lifetimes->at(i);
            } else if (hb_lifetimes->at(i)<min_lifetime) {
                min_lifetime = hb_lifetimes->at(i);
            }
            sum+=hb_lifetimes->at(i);
        }
        double avg_lifetime=1.0*sum/hb_lifetimes->size();

        printf("%d unique hbonds found\n",(hbonds)->size());
        printf("average lifetime: %f\n",avg_lifetime);
        printf("hbonds over time (%d steps):\n", num_timesteps);
        for (int i=0; i<num_timesteps; ++i) {
            printf("step %d, %d hbonds\n",i,num_hbonds[i]);
        }
        printf("Ch pre-normalization:\n");
        for (int i=0; i<max_meas_time; ++i) {
            printf("%f ",Ch[i]);
        }
        printf("\n");
        // normalizing R, Cv, and MSD
        int t0_max;
        //double omega;

        if (VERBOSE) {
            printf("Theta final: [");
            for (int i=0; i<max_meas_time; ++i) {
                int my_2d_index=get_idx(i,1394,num_particles);
                printf("%d  %d\n",my_2d_index,theta[my_2d_index]);
            }
        }

        printf("N_occu_total = %d\n",N_occu_total);
        double N_occu_ave = 1.0*N_occu_total/num_timesteps;
        //omega = (1.0*omega_total/N_occu_ave)/(num_timesteps);
        printf("N_occu_ave = %f\n",N_occu_ave);
        for (int i=0; i<max_meas_time; ++i) {
            t0_max=num_timesteps-i;
            Ch[i]=Ch[i]/t0_max/num_particles; // normalizing Ch
            if (VERBOSE) {
                printf("%d\n",t0_max);
            }

            plot_time[i] = (i*dump_time)/1000.0; // now in ps
        }
        printf("Ch post-normalization:\n");
        for (int i=0; i<max_meas_time; ++i) {
            printf("%f ",Ch[i]);
        }
        printf("\n");
        for (int i=0; i<num_timesteps; ++i) {
            full_time[i] = (i*dump_time)/1000.0; // now in ps
        }

        // FIXME
           ofstream CHout;
           CHout.open("results/C_H-"+suffix+"-w"+to_string(num_particles)+"_s"+to_string(skip)+
                   "_m"+to_string(max_meas_time)+"_t"+to_string(num_timesteps)+"_d"+to_string(dump_time)+"fs");
           ofstream HBout;
           HBout.open("results/H_B-"+suffix+"-w"+to_string(num_particles)+"_s"+to_string(skip)+
                   "_m"+to_string(max_meas_time)+"_t"+to_string(num_timesteps)+"_d"+to_string(dump_time)+"fs");
           CHout << "#t Ch" << endl;
           for (int i=0; i< max_meas_time; ++i) {
               CHout << plot_time[i] << " " << Ch[i] << endl;
           }
           HBout << "#t hbonds" << endl;
           for (int i=0; i< num_timesteps; ++i) {
               HBout << full_time[i] << " " << num_hbonds[i] << endl;
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
    double* Ch;
    double* positions;
    double* h_positions;
    double* plot_time;
    double* full_time;
    double oo_cutoff=5.00; // FIXME read in from arguments
    double oh_cutoff=3.6;
    double phi_cutoff=150;
    double x0=-6.05;
    double x1=6.05;
    double y0=-5.25;
    double y1=5.25;
    double z0=16.1;
    double z1=22.90;
    // FIXME: read in box width and box height from the file
    //double box_width=18.6269*2;
    double box_width=9.29072*2;
    //double box_height=18.6269*2;
    double box_height=9.29072*2;
    int skip_steps=0;
    int num_particles=1000;
    int num_timesteps=200;
    int max_meas_time=100;
    int dump_time=20;
    int O_id=2;
    clock_t start, time;
    string suffix;
    bool VERBOSE;

    // FIXME
    const int mb = ceil(box_width / (BIN_SIZE)); // number of bins per side
    const int nb = ceil(box_height / (BIN_SIZE));
    vector<int>** bins;
    bins = new vector<int>*[(nb * mb * mb)];
    printf("bin size is %f, number of bins: %d\n",BIN_SIZE, nb*mb*mb);
    printf("mb: %d, nb: %d\n",mb,nb);
    for (int i=0; i<mb; ++i) {
        for (int j=0; j<mb; ++j) {
            for (int k=0; k<nb; ++k) {
                tuple<int, int, int> block_id;
                get<0>(block_id) = i;
                get<1>(block_id) = j;
                get<2>(block_id) = k;
                int id = get_bin_idx(block_id, mb, nb);
                bins[id] = new vector<int>();
            }
        }
    }

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
        if (strcmp("-r", argv[i])==0) {
            dump_time = atoi(argv[i+1]);
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
            cout << "Verbose output"<< endl;
        }
        else { VERBOSE=0; }
    }


    //pair<int, int> test(4,5);
    //pair<int, int> testt(6,5);
    //int* testvector;
    //int* testvector2;
    hbond_map hbonds;
    //*hbonds ={
    //    {test,testvector},
    //    {testt,testvector2} };
    //printf("empty? %d",hbonds->empty());
    start=clock();

    read_hbonds(file_name,theta,Ch,
            positions,h_positions,
            bins, mb, nb,
            oo_cutoff, oh_cutoff, phi_cutoff,&hbonds,
            plot_time,full_time,skip_steps,dump_time,
            num_particles,num_timesteps,max_meas_time,
            x0,x1,y0,y1,z0,z1,suffix,O_id,VERBOSE);
    time=clock()-start;
    time=time/CLOCKS_PER_SEC;
    cout << "Took " << time << " sec to run" << endl;
    cout << "--------------------" << endl;
    return 0;
}
