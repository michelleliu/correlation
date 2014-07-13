#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"

#include <functional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#define BIN_SIZE 5*cutoff

using namespace std;

struct pairhash{
  size_t operator()(const pair<int, int> &p) const {
      return
        hash<unsigned int>()(p.first) ^
        hash<unsigned int>()(p.second);
  }
};


struct clump_t {
  vector<particle_t>* particles;
  vector<int>** bins;
  vector<int>** ghost_bins;
  omp_lock_t* locks;
  int clump_id;
  int num_bins;
  double clump_startx;
  double clump_starty;
};

int get_ghost_idx(int row, int col, const int num) {
  return (row + 1) * (num + 2) + (col + 1);
}

pair<int, int> get_ghost_row_col(int idx, const int num) {
  int col = idx % (num + 2) - 1;
  int row = idx / (num + 2) - 1;
  return make_pair(row, col);
}

bool is_this_ghost_real(){
  //I believe in ghosts!!
  return true;
}

int get_idx(int row, int col, const int num) {
  return (row) * (num) + (col);
}

pair<int, int> get_row_col(int idx, const int num) {
  int col = idx % (num) ;
  int row = idx / (num) ;
  return make_pair(row, col);
}

pair<int, double> get_num_clumps(double size, int nprocs) {
  int num_clumps = (int)(ceil(sqrt((double)nprocs)));
  double clump_width = size / num_clumps;
  return make_pair(num_clumps, clump_width);
}

int get_clump(particle_t particle, double clump_width, int num_clumps) {
  int row_block = (int) floor(particle.x / clump_width);
  int col_block = (int) floor(particle.y / clump_width);
  return get_idx(row_block, col_block, num_clumps);
}

bool are_you_a_ghost(pair<int, int> block_id, int num_bins) {
  return !(block_id.first >= 0 && block_id.first < num_bins &&
          block_id.second >= 0 && block_id.second < num_bins);
}

clump_t* get_clump_from_id(vector<clump_t*>* clumps, int clump_id) {
  for (int i = 0; i < clumps->size(); ++i) {
    if (clumps->at(i)->clump_id == clump_id) {
      return clumps->at(i);
    }
  }
  return NULL;
}


void clear_bins(vector<int>** bins, omp_lock_t* locks, const int num_bins) {
  for (int i = 0; i < (num_bins + 2)*(num_bins+2); ++i) {
    //omp_set_lock(&locks[i]);
    bins[i]->clear();
    //omp_unset_lock(&locks[i]);
  }
}

void populate_bins(vector<int>** bins,
    vector<particle_t>* particles, omp_lock_t* locks, int num_bins,
    double clump_startx, double clump_starty, int startIdx, double clump_width) {
  // printf("Populating %d particles\n", particles->size() - startIdx);
  // #pragma omp for
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  for (int i = startIdx; i < particles->size(); ++i) {
    // Get x, y co-ordinates of this particle
    int row_block = (int) floor((particles->at(i).x - clump_startx) / (BIN_SIZE));
    int col_block = (int) floor((particles->at(i).y - clump_starty) / (BIN_SIZE));

    // Force this particle into a ghost bin if its beyond our clump
    if (particles->at(i).x >= (clump_startx + clump_width) && row_block != num_bins) {
      row_block = num_bins;
    }

    if (particles->at(i).y >= (clump_starty + clump_width) && col_block != num_bins) {
      col_block = num_bins;
    }
    // if (startIdx != 0) {
    //   printf("DBG Rank: %d x %lf y %lf startx %lf starty %lf, row %d, col %d\n",
    //     rank, particles->at(i).x, particles->at(i).y, clump_startx, clump_starty, row_block, col_block);
    // }

    int id = get_ghost_idx(row_block, col_block, num_bins);
    omp_set_lock(&locks[id]);
    bins[id]->push_back(i);
    omp_unset_lock(&locks[id]);
  }
}

//place particles that belong to a specific clumps into bins
void bin_clump(clump_t* clump, double clump_width){
  //initialize the bins
  int num_bins = ceil(clump_width/(BIN_SIZE));
  clump->num_bins = num_bins;
  clump->bins = new vector<int>*[(num_bins+2)*(num_bins+2)];
  clump->locks = new omp_lock_t[(num_bins+2)*(num_bins+2)];
  for(int i=-1; i< num_bins+1; ++i){
    for (int j=-1; j < num_bins+1; ++j){
      int id=get_ghost_idx(i,j,num_bins);
      clump->bins[id] = new vector<int>();
      omp_init_lock(&clump->locks[id]);
    }
  }
  //populate bin
  populate_bins(clump->bins, clump->particles, clump->locks, num_bins,
                clump->clump_startx, clump->clump_starty, 0, clump_width);
}

void get_neighbors(const pair<int, int>& block_id, vector<pair<int, int> >* neighbors) {
  // printf("For block (%d, %d) nbrs: ", block_id.first, block_id.second);
  for (int i = block_id.first - 1; i <= block_id.first + 1; ++i) {
    for (int j = block_id.second - 1; j <= block_id.second + 1; ++j) {
      if (!(i == block_id.first && j == block_id.second)) {
        // printf(" (%d, %d), ", i, j);
        neighbors->push_back(make_pair(i, j));
      }
    }
  }
}

void get_neighbor_clumps(const pair<int, int>& clump_id, vector<pair<int, int> >* neighbors, int num_clumps) {
  for (int i = clump_id.first - 1; i <= clump_id.first + 1; ++i) {
    for (int j = clump_id.second - 1; j <= clump_id.second + 1; ++j) {
      if (i >= 0 && j >= 0 && i < num_clumps && j < num_clumps) {
        if (!(i == clump_id.first && j == clump_id.second)) {
          // printf(" (%d, %d), ", i, j);
          neighbors->push_back(make_pair(i, j));
        }
      }
    }
  }
}

void get_particles_from_bins(vector<pair<int, int> >* bins_to_copy,
    vector<int>** my_bins,
    vector<particle_t>* my_particles,
    int num_bins,
    vector<particle_t>* all_particles) {
  for (int i = 0; i < bins_to_copy->size(); ++i) {
    vector<int>* bin = my_bins[get_ghost_idx(bins_to_copy->at(i).first,
        bins_to_copy->at(i).second, num_bins)];
    for (int j = 0; j < bin->size(); ++j) {
      all_particles->push_back(my_particles->at(bin->at(j)));
    }
  }
}

void add_particles_to_clump_arr(particle_t* particles_to_add, int size, 
    vector<particle_t>* my_particles) {
  for (int i = 0; i < size; ++i) {
    my_particles->push_back(particles_to_add[i]);
  }
}

void add_particles_to_clump(vector<particle_t>& particles_to_add, vector<particle_t>* my_particles) {
  for (int i = 0; i < particles_to_add.size(); ++i) {
    my_particles->push_back(particles_to_add.at(i));
  }
}

void get_ghost_bins(const pair<int, int>& my_clump_id, const pair<int, int>& nbr_clump_id,
    vector<pair<int, int> >* send_bins, vector<pair<int, int> >* recv_bins, int num_bins) {
  int clump_diff = 10*(my_clump_id.first - nbr_clump_id.first) + (my_clump_id.second - nbr_clump_id.second);
  switch (clump_diff) {
    case 11:
      send_bins->push_back(make_pair(0,0));
      recv_bins->push_back(make_pair(-1,-1));
      break;
    case 9:
      send_bins->push_back(make_pair(0,num_bins-1));
      recv_bins->push_back(make_pair(-1,num_bins));
      break;
    case -11:
      send_bins->push_back(make_pair(num_bins-1,num_bins-1));
      recv_bins->push_back(make_pair(num_bins,num_bins));
      break;
    case -9:
      send_bins->push_back(make_pair(num_bins-1,0));
      recv_bins->push_back(make_pair(num_bins,-1));
      break;
    case 10:
      for (int i = 0; i < num_bins; ++i) {
        send_bins->push_back(make_pair(0,i));
        recv_bins->push_back(make_pair(-1,i));
      }
      break;
    case -10:
      for (int i = 0; i < num_bins; ++i) {
        send_bins->push_back(make_pair(num_bins-1,i));
        recv_bins->push_back(make_pair(num_bins,i));
      }
      break;
    case 1:
      for (int i = 0; i < num_bins; ++i) {
        send_bins->push_back(make_pair(i,0));
        recv_bins->push_back(make_pair(i,-1));
      }
      break;
    case -1:
      for (int i = 0; i < num_bins; ++i) {
        send_bins->push_back(make_pair(i,num_bins-1));
        recv_bins->push_back(make_pair(i,num_bins));
      }
      break;
    default:
      printf("GHOSTS");
  }
}

void apply_force_bins(vector<int>* source, vector<int>* other, vector<particle_t>* particles,
    double *dmin, double *davg, int *navg) {
  for (int iaf = 0; iaf < source->size(); ++iaf) {
    for (int jaf = 0; jaf < other->size(); ++jaf) {
      // printf("Applying force to %d from %d\n", source->at(i), other->at(j));
      apply_force(particles->at(source->at(iaf)),
                  particles->at(other->at(jaf)),
                  dmin, davg, navg);
    }
  }
}

int get_pid(int clump_id, int nproc)
{
  return clump_id % nproc; 
}

bool is_particle_in_clump(particle_t p, clump_t* clump, double clump_width) {
  // Particle is in clump if its position is within clump boundaries
  if (p.x >= clump->clump_startx && p.x < clump->clump_startx + clump_width &&
      p.y >= clump->clump_starty && p.y < clump->clump_starty + clump_width) {
    return true;
  }
  return false;
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{
  int navg, nabsavg=0;
  double dmin, absmin=1.0,davg,absavg=0.0;
  double rdavg,rdmin;
  int rnavg;

  //
  //  process command line parameters
  //
  if( find_option( argc, argv, "-h" ) >= 0 )
  {
    printf( "Options:\n" );
    printf( "-h to see this help\n" );
    printf( "-n <int> to set the number of particles\n" );
    printf( "-o <filename> to specify the output file name\n" );
    printf( "-s <filename> to specify a summary file name\n" );
    printf( "-no turns off all correctness checks and particle output\n");
    return 0;
  }

  int n = read_int( argc, argv, "-n", 1000 );
  char *savename = read_string( argc, argv, "-o", NULL );
  char *sumname = read_string( argc, argv, "-s", NULL );

  //
  //  set up MPI
  //
  int n_proc, rank;
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  //
  //  allocate generic resources
  //
  FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
  FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;

  MPI_Datatype PARTICLE;
  MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
  MPI_Type_commit( &PARTICLE );

  //
  //  initialize and distribute the particles (that's fine to leave it unoptimized)
  //
  set_size( n );
  particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

  if (rank == 0) {
    init_particles( n, particles );
  }

  pair<int, double> clump_info = get_num_clumps(get_size(), n_proc);
  int num_clumps = clump_info.first;
  double clump_width = clump_info.second;
  const int num_bins = ceil(clump_width / (BIN_SIZE));
  if (rank == 0) {
    printf("num_clumps %d with width %lf\n", num_clumps, clump_width);
    printf("Size is %lf, bin size is %lf, Number of bins %d\n", get_size(), BIN_SIZE,
        num_bins);
  }

  vector<clump_t*> my_clumps;

  // Assign clumps to processors
  for (int i = 0; i < (num_clumps) * (num_clumps); ++i) {
    if ( i % n_proc == rank ) {
      clump_t* c = new clump_t;
      c->clump_id = i;
      pair<int, int> clump_id_pair = get_row_col(i, num_clumps);
      c->clump_startx = clump_id_pair.first * clump_width;
      c->clump_starty = clump_id_pair.second * clump_width;
      c->particles = new vector<particle_t>;
      my_clumps.push_back(c);
    }
  }

  MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);
  //Assign particles to clumps
  for (int i = 0; i < n; ++i) {
    int clump_id = get_clump(particles[i], clump_width, num_clumps);
    if (clump_id % n_proc == rank) {
      for (int j=0; j< my_clumps.size(); j++){
        if (my_clumps[j]->clump_id == clump_id){
          my_clumps[j]->particles->push_back(particles[i]);
        }
      }
    }
  }
/*------------------ Done distributing particles ----------------*/

  //printf("DONE DISTRIBUTING PARTICLES!\n");

  //
  //  simulate a number of time steps
  //
  double simulation_time = read_timer( );
  #pragma omp parallel private(dmin)
  {
  //Do binning in each processor
  for(int i=0; i<my_clumps.size(); ++i){
    bin_clump(my_clumps[i], clump_width);
  }

  for( int step = 0; step < NSTEPS; step++ )
  {
    navg = 0;
    dmin = 1.0;
    davg = 0.0;
    // printf("Start step %d\n", step);

    //
    //  save current step if necessary (slightly different semantics than in other codes)
    //
    if( find_option( argc, argv, "-no" ) == -1 )
      if( fsave && (step%SAVEFREQ) == 0 )
        save( fsave, n, particles );

    //
    //  Exchange ghost bins between processors
    //
    unordered_map<pair<int, int>, pair<particle_t*, int>, pairhash> self_clump_transfers;
  
    // Map from clump index (i.e. k) to various datastructures. 
    unordered_map<int, particle_t**> particles_to_send_map;
    unordered_map<int, particle_t**> particles_to_recv_map;

    unordered_map<int, int*> send_sizes_map;
    unordered_map<int, int*> recv_sizes_map;
    unordered_map<int, MPI_Request*> size_reqs_map;
    unordered_map<int, MPI_Request*> reqs_map;

    unordered_map<int, MPI_Status*> size_status_map;
    unordered_map<int, MPI_Status*> status_map;

    for (int k=0; k < my_clumps.size(); ++k) {
      // printf("BEGIN Rank: %d clump: %i has %d particles\n", rank, k, my_clumps[k]->particles->size());
      vector<pair<int,int> > nbr_clumps;
      pair<int, int> my_clump_id = get_row_col(my_clumps[k]->clump_id, num_clumps);
      get_neighbor_clumps(my_clump_id, &nbr_clumps, num_clumps);

      vector<pair<int, int> > send_bins;
      vector<pair<int, int> > recv_bins;

      int num_bins = my_clumps[k]->num_bins;
      vector<int>** bins = my_clumps[k]->bins;
      omp_lock_t* locks = my_clumps[k]->locks;
      vector<particle_t>* particles = my_clumps[k]->particles;

      MPI_Request* reqs = new MPI_Request[2*nbr_clumps.size()];
      reqs_map[k] = reqs;
      MPI_Request* size_reqs = new MPI_Request[2*nbr_clumps.size()];
      size_reqs_map[k] = size_reqs;

      MPI_Status* status = new MPI_Status[2*nbr_clumps.size()];
      status_map[k] = status;
      MPI_Status* size_status  = new MPI_Status[2*nbr_clumps.size()];
      size_status_map[k] = size_status;

      particle_t** particles_to_send = new particle_t*[nbr_clumps.size()];
      particles_to_send_map[k] = particles_to_send;

      int* send_sizes = new int[nbr_clumps.size()];
      send_sizes_map[k] = send_sizes;
      int* recv_sizes = new int[nbr_clumps.size()];
      recv_sizes_map[k] = recv_sizes;

      // printf("FINISHED MPI SETUP\n");

      vector<particle_t> particles_temp;

      int my_proc = get_pid(my_clumps[k]->clump_id, n_proc);
      for (int j=0; j < nbr_clumps.size(); ++j) {
        get_ghost_bins(my_clump_id, nbr_clumps[j], &send_bins, &recv_bins, num_bins);
        int nbr_clump_id = get_idx(nbr_clumps[j].first,nbr_clumps[j].second, num_clumps);
        int other_proc = get_pid(nbr_clump_id, n_proc);

        // particles_to_send[j] = new vector<particle_t>();
        get_particles_from_bins(&send_bins, bins,
            particles, num_bins, &particles_temp);

        send_sizes[j] = particles_temp.size();
        // printf("Rank: %d, Sending size %d from c %d to c %d\n", rank, send_sizes[j],
        //         my_clumps[k]->clump_id, nbr_clump_id);
        particles_to_send[j] = new particle_t[send_sizes[j]];
        copy(particles_temp.begin(), particles_temp.end(), particles_to_send[j]);

        if (my_proc < other_proc) {
          MPI_Isend(&send_sizes[j], 1, MPI_INT, other_proc, 0, MPI_COMM_WORLD, &size_reqs[2*j]);
          MPI_Irecv(&recv_sizes[j], 1, MPI_INT, other_proc, 0, MPI_COMM_WORLD, &size_reqs[2*j+1]);
        }
        else if (my_proc > other_proc) {
          MPI_Irecv(&recv_sizes[j], 1, MPI_INT, other_proc, 0, MPI_COMM_WORLD, &size_reqs[2*j+1]);
          MPI_Isend(&send_sizes[j], 1, MPI_INT, other_proc, 0, MPI_COMM_WORLD, &size_reqs[2*j]);
        } else {
          size_reqs[2*j+1] = MPI_REQUEST_NULL;
          size_reqs[2*j] = MPI_REQUEST_NULL;
          // printf("ADDING Self transfer to %d from %d 2*j %d 2*j+1 %d\n", nbr_clump_id, my_clumps[k]->clump_id, 2*j, 2*j+1);
          self_clump_transfers.insert(make_pair(
                make_pair(nbr_clump_id, my_clumps[k]->clump_id),
                make_pair(particles_to_send[j], send_sizes[j])));
        }
        send_bins.clear();
        recv_bins.clear();
        particles_temp.clear();
      }
    }

    for (int k = 0; k < my_clumps.size(); ++k) {
      vector<pair<int,int> > nbr_clumps;
      pair<int, int> my_clump_id = get_row_col(my_clumps[k]->clump_id, num_clumps);
      get_neighbor_clumps(my_clump_id, &nbr_clumps, num_clumps);

      MPI_Waitall(2*nbr_clumps.size(), size_reqs_map[k], size_status_map[k]);
    }

    // printf("FINISHED SIZE EXCHANGE\n");


    for (int k=0; k < my_clumps.size(); ++k) {
      vector<pair<int,int> > nbr_clumps;
      pair<int, int> my_clump_id = get_row_col(my_clumps[k]->clump_id, num_clumps);
      get_neighbor_clumps(my_clump_id, &nbr_clumps, num_clumps);

      vector<pair<int, int> > send_bins;
      vector<pair<int, int> > recv_bins;

      int num_bins = my_clumps[k]->num_bins;
      vector<int>** bins = my_clumps[k]->bins;
      omp_lock_t* locks = my_clumps[k]->locks;
      vector<particle_t>* particles = my_clumps[k]->particles;

      particle_t** particles_to_recv = new particle_t*[nbr_clumps.size()];
      particles_to_recv_map[k] = particles_to_recv;
      int* send_sizes = send_sizes_map[k];
      int* recv_sizes = recv_sizes_map[k];
      particle_t** particles_to_send = particles_to_send_map[k];
      MPI_Request* reqs = reqs_map[k];

      int my_proc = get_pid(my_clumps[k]->clump_id, n_proc);
    
      for (int j=0; j < nbr_clumps.size(); ++j) {
        int nbr_clump_id = get_idx(nbr_clumps[j].first,nbr_clumps[j].second, num_clumps);
        int other_proc = get_pid(nbr_clump_id, n_proc);

        // particles_to_recv[j] = new vector<particle_t>(recv_sizes[j]);
        if (my_proc < other_proc) {
          particles_to_recv[j] = new particle_t[recv_sizes[j]];

          // printf("Rank: %d, j %d, nbr %d got size from c %d at p %d of %d\n", rank, j, nbr_clump_id, my_clumps[k]->clump_id, other_proc, recv_sizes[j]); 
          // printf("Rank: %d, j %d, Sending size %d\n", rank, j, send_sizes[j]);
          MPI_Isend(particles_to_send[j], send_sizes[j],
              PARTICLE, other_proc, 0, MPI_COMM_WORLD,  &reqs[2*j]);
          MPI_Irecv(particles_to_recv[j], recv_sizes[j],
              PARTICLE, other_proc, 0, MPI_COMM_WORLD,  &reqs[2*j+1]);
        }
        else if (my_proc > other_proc) {
          particles_to_recv[j] = new particle_t[recv_sizes[j]];
          // printf("Rank: %d, j %d, nbr %d got size from c %d at p %d of %d\n", rank, j, nbr_clump_id, my_clumps[k]->clump_id, other_proc, recv_sizes[j]); 
          // printf("Rank: %d, j %d, Sending size %d\n", rank, j, send_sizes[j]);
          MPI_Irecv(particles_to_recv[j], recv_sizes[j],
              PARTICLE, other_proc, 0, MPI_COMM_WORLD,  &reqs[2*j+1]);
          MPI_Isend(particles_to_send[j], send_sizes[j],
              PARTICLE, other_proc, 0, MPI_COMM_WORLD,  &reqs[2*j]);
        } else {
          reqs[2*j+1] = MPI_REQUEST_NULL;
          reqs[2*j] = MPI_REQUEST_NULL;
          // if (self_clump_transfers.find(make_pair(my_clumps[k]->clump_id, nbr_clump_id)) == self_clump_transfers.end()) {
          //   printf("ERROR Self transfer to %d from %d not found\n", my_clumps[k]->clump_id, nbr_clump_id);
          // } 
          // pair<particle_t*, int> recv =
          //   self_clump_transfers.find(make_pair(my_clumps[k]->clump_id, nbr_clump_id))->second;
          // printf("DOING A SELF TRANSFER of size %d\n", recv.second);
          // copy(recv.first, recv.first + recv.second, particles_to_recv[j]);
          // recv_sizes[j] = recv.second;
          recv_sizes[j] = 0;
        }
      }
      // printf("AFTER  Rank: %d clump: %i has %d particles\n", rank, k, my_clumps[k]->particles->size());
    }

    for (int k=0; k < my_clumps.size(); ++k) {
      vector<pair<int,int> > nbr_clumps;
      pair<int, int> my_clump_id = get_row_col(my_clumps[k]->clump_id, num_clumps);
      get_neighbor_clumps(my_clump_id, &nbr_clumps, num_clumps);

      int num_bins = my_clumps[k]->num_bins;
      vector<int>** bins = my_clumps[k]->bins;
      omp_lock_t* locks = my_clumps[k]->locks;
      vector<particle_t>* particles = my_clumps[k]->particles;

      int* recv_sizes = recv_sizes_map[k];
      particle_t** particles_to_recv = particles_to_recv_map[k];

      MPI_Waitall(2*nbr_clumps.size(), reqs_map[k], status_map[k]);
      // printf("FINISHED PARTICLE EXCHANGE\n");
      // MPI_Barrier(MPI_COMM_WORLD);

      for (int j = 0; j < nbr_clumps.size(); ++j) {
        int oldSize = particles->size();
        if (recv_sizes[j] > 0) {
          // printf("Rank: %d j: %d size of ptrecv %d\n", rank, j, particles_to_recv[j]->size());
          add_particles_to_clump_arr(particles_to_recv[j], recv_sizes[j], particles);

          populate_bins(bins,
              particles, locks, num_bins,
              my_clumps[k]->clump_startx, my_clumps[k]->clump_starty, oldSize, clump_width);
          delete particles_to_recv[j];
        }
      }
    }

    // DO the self transfers here
    for(unordered_map<pair<int, int>, pair<particle_t*, int>, pairhash>::iterator itr = self_clump_transfers.begin(); itr != self_clump_transfers.end(); ++itr) {
      clump_t* dest_clump = get_clump_from_id(&my_clumps, itr->first.first);
      vector<particle_t>* particles = dest_clump->particles;
      vector<int>** bins = dest_clump->bins;
      omp_lock_t* locks = dest_clump->locks;
      int num_bins = dest_clump->num_bins;

      int oldSize = particles->size();
      add_particles_to_clump_arr(itr->second.first, itr->second.second, particles);

      populate_bins(bins,
          particles, locks, num_bins,
          dest_clump->clump_startx, dest_clump->clump_starty, oldSize, clump_width);
      
      delete itr->second.first;
    }

    self_clump_transfers.clear();

    for (int k = 0; k < my_clumps.size(); ++k) {
        delete[] particles_to_send_map[k];
        delete[] particles_to_recv_map[k];
        delete[] send_sizes_map[k];
        delete[] recv_sizes_map[k];
        delete[] size_reqs_map[k];
        delete[] reqs_map[k];
        delete[] size_status_map[k];
        delete[] status_map[k];
    }

    // printf("FINISHED MPI EXCHANGE ?\n");

    // compute all forces
    // for bin in bins
    //    get neigboring bins
    //    for each neigboring bin
    //      do apply_force
    #pragma omp for reduction (+:navg) reduction(+:davg) collapse(2)
    for (int k=0; k< my_clumps.size(); ++k){
      for (int j = 0; j< (num_bins+2)*(num_bins+2); ++j) {
        int num_bins = my_clumps[k]->num_bins;
        vector<int>** bins = my_clumps[k]->bins;
        vector<particle_t>* particles = my_clumps[k]->particles;

        vector<int>* this_bin = bins[j];
        // List surrounding bins using key
        pair<int, int> block_id = get_ghost_row_col(j, num_bins);

        if (!are_you_a_ghost(block_id, num_bins)) {
          for (int i = 0; i < this_bin->size(); ++i) {
            particles->at(this_bin->at(i)).ax = particles->at(this_bin->at(i)).ay = 0;
          }

          apply_force_bins(this_bin, this_bin, particles, &dmin, &davg, &navg);
          vector<pair<int, int> > neighbors;
          get_neighbors(block_id, &neighbors);

          for (vector<pair<int, int> >::iterator nbr_itr = neighbors.begin();
               nbr_itr != neighbors.end(); ++nbr_itr) {
            vector<int>* nbr_bin = bins[get_ghost_idx(nbr_itr->first, nbr_itr->second, num_bins)];
            apply_force_bins(this_bin, nbr_bin, particles, &dmin, &davg, &navg);
          }
        }
      }
    }

    // printf("FINISHED APPLY FORCE !!\n");

    if( find_option( argc, argv, "-no" ) == -1 )
    {

      MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);


      if (rank == 0){
        //
        // Computing statistical data
        //
        if (rnavg) {
          absavg +=  rdavg/rnavg;
          nabsavg++;
        }
        if (rdmin < absmin) absmin = rdmin;
      }
    }

    //
    //  move particles
    //
    for( int i = 0; i < my_clumps.size(); i++ ){
      int num_bins = my_clumps[i]->num_bins;
      vector<int>** bins = my_clumps[i]->bins;

      #pragma omp for
      for (int j = 0; j< (num_bins+2)*(num_bins+2); ++j) {
        vector<int>* this_bin = bins[j];
        // List surrounding bins using key
        pair<int, int> block_id = get_ghost_row_col(j, num_bins);

        if (!are_you_a_ghost(block_id, num_bins)) {
          for (int k = 0; k < this_bin->size(); ++k) {
            move( my_clumps[i]->particles->at(this_bin->at(k)));
          }
        }
      }
    }

    // printf("FINISHED MOVE !!\n");

    // Figure out which particles moved out of our clump and which came in ?
    //
    // Exchange this with those processors
    unordered_map<pair<int, int>, vector<particle_t>, pairhash > particles_to_move;
    for (int i = 0; i < my_clumps.size(); ++i) {
      int num_bins = my_clumps[i]->num_bins;
      int my_clump_id = my_clumps[i]->clump_id;
      vector<int>** bins = my_clumps[i]->bins;
      vector<particle_t> new_particles_in_clump;
      for (int j = 0; j< (num_bins+2)*(num_bins+2); ++j) {
        vector<int>* this_bin = bins[j];
        pair<int, int> block_id = get_ghost_row_col(j, num_bins);

        if (!are_you_a_ghost(block_id, num_bins)) {
          for (int k = 0; k < this_bin->size(); ++k) {
            if (is_particle_in_clump(my_clumps[i]->particles->at(this_bin->at(k)), my_clumps[i], clump_width)) {
              new_particles_in_clump.push_back(my_clumps[i]->particles->at(this_bin->at(k)));
            } else {
              int dest_clump_id = get_clump(my_clumps[i]->particles->at(this_bin->at(k)),
                  clump_width, num_clumps);
              vector<particle_t>& parts = particles_to_move[make_pair(dest_clump_id, my_clump_id)];
              parts.push_back(my_clumps[i]->particles->at(this_bin->at(k)));
            }
          }
        }
      }
      // Clear all our old particles and put new particles in it.
      // printf("PRE-SWAP Rank: %d clump: %i has %d particles\n", rank, i, my_clumps[i]->particles->size());
      my_clumps[i]->particles->clear();
      my_clumps[i]->particles->swap(new_particles_in_clump);
      // printf("POST-SWAP Rank: %d clump: %i has %d particles\n", rank, i, my_clumps[i]->particles->size());
    }

    //for (unordered_map<pair<int, int>, vector<particle_t>, pairhash>::iterator itr = particles_to_move.begin(); itr != particles_to_move.end(); ++itr) {
      // printf("Rank: %d, clump %d sending %d to %d\n", rank, itr->first.second, itr->second.size(), itr->first.first);
    //}


    particles_to_send_map.clear();
    particles_to_recv_map.clear();

    send_sizes_map.clear();
    recv_sizes_map.clear();
    size_reqs_map.clear();
    reqs_map.clear();

    size_status_map.clear();
    status_map.clear();

    // Send, receive particles to neighbor clumps
    for (int i = 0; i < my_clumps.size(); ++i) {
      vector<particle_t>* particles = my_clumps[i]->particles;

      // Get the neigbor clumps
      vector<pair<int,int> > nbr_clumps;
      pair<int, int> my_clump_id = get_row_col(my_clumps[i]->clump_id, num_clumps);
      get_neighbor_clumps(my_clump_id, &nbr_clumps, num_clumps);

      MPI_Request* reqs = new MPI_Request[2*nbr_clumps.size()];
      reqs_map[i] = reqs;

      MPI_Request* size_reqs = new MPI_Request[2*nbr_clumps.size()];
      size_reqs_map[i] = size_reqs;

      MPI_Status* status = new MPI_Status[2*nbr_clumps.size()];
      status_map[i] = status;
      MPI_Status* size_status  = new MPI_Status[2*nbr_clumps.size()];
      size_status_map[i] = size_status;

      int* send_sizes = new int[nbr_clumps.size()];
      send_sizes_map[i] = send_sizes;
      int* recv_sizes = new int[nbr_clumps.size()];
      recv_sizes_map[i] = recv_sizes;

      int my_proc = get_pid(my_clumps[i]->clump_id, n_proc);

      for (int j=0; j < nbr_clumps.size(); ++j) {
        int nbr_clump_id = get_idx(nbr_clumps[j].first,nbr_clumps[j].second, num_clumps);
        int other_proc = get_pid(nbr_clump_id, n_proc);

        send_sizes[j] = 0;
        if (particles_to_move.find(make_pair(nbr_clump_id, my_clumps[i]->clump_id)) != 
            particles_to_move.end()) {
          send_sizes[j] = particles_to_move.find(
            make_pair(nbr_clump_id, my_clumps[i]->clump_id))->second.size();
        }

        if (my_proc < other_proc) {
          MPI_Isend(&send_sizes[j], 1, MPI_INT, other_proc, 0, MPI_COMM_WORLD, &size_reqs[2*j]);
          MPI_Irecv(&recv_sizes[j], 1, MPI_INT, other_proc, 0, MPI_COMM_WORLD, &size_reqs[2*j+1]);
        }
        else if (my_proc > other_proc) {
          MPI_Irecv(&recv_sizes[j], 1, MPI_INT, other_proc, 0, MPI_COMM_WORLD, &size_reqs[2*j+1]);
          MPI_Isend(&send_sizes[j], 1, MPI_INT, other_proc, 0, MPI_COMM_WORLD, &size_reqs[2*j]);
        } else {
          size_reqs[2*j] = MPI_REQUEST_NULL;
          size_reqs[2*j+1] = MPI_REQUEST_NULL;
        }
      }
    }

    for (int k = 0; k < my_clumps.size(); ++k) {
      vector<pair<int,int> > nbr_clumps;
      pair<int, int> my_clump_id = get_row_col(my_clumps[k]->clump_id, num_clumps);
      get_neighbor_clumps(my_clump_id, &nbr_clumps, num_clumps);

      MPI_Waitall(2*nbr_clumps.size(), size_reqs_map[k], size_status_map[k]);
    }

    // printf("FINISHED TRANSFERRING SIZES\n");

    for (int k=0; k < my_clumps.size(); ++k) {
      vector<pair<int,int> > nbr_clumps;
      pair<int, int> my_clump_id = get_row_col(my_clumps[k]->clump_id, num_clumps);
      get_neighbor_clumps(my_clump_id, &nbr_clumps, num_clumps);

      int num_bins = my_clumps[k]->num_bins;
      vector<int>** bins = my_clumps[k]->bins;
      omp_lock_t* locks = my_clumps[k]->locks;
      vector<particle_t>* particles = my_clumps[k]->particles;

      particle_t** particles_to_recv = new particle_t*[nbr_clumps.size()];
      particles_to_recv_map[k] = particles_to_recv;
      int* send_sizes = send_sizes_map[k];
      int* recv_sizes = recv_sizes_map[k];
      MPI_Request* reqs = reqs_map[k];
      int my_proc = get_pid(my_clumps[k]->clump_id, n_proc);

      for (int j=0; j < nbr_clumps.size(); ++j) {
        int nbr_clump_id = get_idx(nbr_clumps[j].first,nbr_clumps[j].second, num_clumps);
        int other_proc = get_pid(nbr_clump_id, n_proc);

        particle_t* particles_to_send = NULL;

        if (particles_to_move.find(make_pair(nbr_clump_id, my_clumps[k]->clump_id)) != 
            particles_to_move.end()) {
          particles_to_send = &(particles_to_move.find(
                make_pair(nbr_clump_id, my_clumps[k]->clump_id))->second)[0];
        }
        // printf("TRANSFER Rank: %d, j %d, got size from %d %d\n", rank, j, other_proc, recv_sizes[j]); 
        // printf("TRANSFER Rank: %d, j %d, Sending size %d\n", rank, j, send_sizes[j]);

        if (my_proc < other_proc) {
          particles_to_recv[j] = new particle_t[recv_sizes[j]];
          //if (send_sizes[j] > 0) {
            MPI_Isend(particles_to_send, send_sizes[j],
                PARTICLE, other_proc, 0, MPI_COMM_WORLD, &reqs[2*j]);
          //}

          //if (recv_sizes[j] > 0) {
            MPI_Irecv(particles_to_recv[j], recv_sizes[j],
                PARTICLE, other_proc, 0, MPI_COMM_WORLD, &reqs[2*j+1]);
          //}
        }
        else if (my_proc > other_proc) {
          particles_to_recv[j] = new particle_t[recv_sizes[j]];
          //if (recv_sizes[j] > 0) {
            MPI_Irecv(particles_to_recv[j], recv_sizes[j],
                PARTICLE, other_proc, 0, MPI_COMM_WORLD, &reqs[2*j+1]);
          //}
          //if (send_sizes[j] > 0) {
            MPI_Isend(particles_to_send, send_sizes[j],
                PARTICLE, other_proc, 0, MPI_COMM_WORLD, &reqs[2*j]);
          //}
        } else {
          reqs[2*j] = MPI_REQUEST_NULL;
          reqs[2*j+1] = MPI_REQUEST_NULL;
        
          if (particles_to_move.find(make_pair(my_clumps[k]->clump_id, nbr_clump_id)) != 
              particles_to_move.end()) {
            vector<particle_t>& myvec = particles_to_move.find(
                make_pair(my_clumps[k]->clump_id, nbr_clump_id))->second;
            particles_to_recv[j] = new particle_t[myvec.size()];
            copy( myvec.begin(), myvec.end(), particles_to_recv[j]);
            recv_sizes[j] = myvec.size();
          } else {
            recv_sizes[j] = 0;
          }
        }
      }
    }

    for (int k=0; k < my_clumps.size(); ++k) {
      vector<pair<int,int> > nbr_clumps;
      pair<int, int> my_clump_id = get_row_col(my_clumps[k]->clump_id, num_clumps);
      get_neighbor_clumps(my_clump_id, &nbr_clumps, num_clumps);

      int num_bins = my_clumps[k]->num_bins;
      vector<int>** bins = my_clumps[k]->bins;
      omp_lock_t* locks = my_clumps[k]->locks;
      vector<particle_t>* particles = my_clumps[k]->particles;

      int* recv_sizes = recv_sizes_map[k];
      particle_t** particles_to_recv = particles_to_recv_map[k];
      MPI_Request* reqs = reqs_map[k];
      MPI_Status* status = status_map[k];

      MPI_Waitall(2*nbr_clumps.size(), reqs, status);
      // MPI_Barrier(MPI_COMM_WORLD);
      // printf("FINISHED TRANSFERRING PARTICLES\n");

      for (int j = 0; j < nbr_clumps.size(); ++j) {
        if (recv_sizes[j] > 0) {
          add_particles_to_clump_arr(particles_to_recv[j], recv_sizes[j], particles);
          delete particles_to_recv[j];
        }
      }
      // printf("POST-TRANSFER Rank: %d clump: %i has %d particles\n", rank, i, my_clumps[i]->particles->size());
    }
    particles_to_move.clear();
    // printf("FINISHED ADDING AS WELL\n");

    for (int k = 0; k < my_clumps.size(); ++k) {
        delete[] particles_to_recv_map[k];
        delete[] send_sizes_map[k];
        delete[] recv_sizes_map[k];
        delete[] size_reqs_map[k];
        delete[] reqs_map[k];
        delete[] size_status_map[k];
        delete[] status_map[k];
    }

    #pragma omp master
    {
      for( int i = 0; i < my_clumps.size(); i++ ){
        int num_bins = my_clumps[i]->num_bins;
        vector<int>** bins = my_clumps[i]->bins;
        omp_lock_t* locks = my_clumps[i]->locks;

        clear_bins(bins, locks, num_bins);
      }
    }
    // printf("FINISHED CLEARING BINS\n");
    #pragma omp barrier
    {
      for( int i = 0; i < my_clumps.size(); i++ ){
        int num_bins = my_clumps[i]->num_bins;
        vector<int>** bins = my_clumps[i]->bins;
        omp_lock_t* locks = my_clumps[i]->locks;
        vector<particle_t>* particles = my_clumps[i]->particles;

        populate_bins(bins, particles, locks, num_bins, my_clumps[i]->clump_startx,
            my_clumps[i]->clump_starty, 0, clump_width);
      }
    }
    // printf("FINISHED RE-BINNING BINS\n");
    // exit(0);
    
    MPI_Barrier(MPI_COMM_WORLD);
  }
  }
  simulation_time = read_timer( ) - simulation_time;

  if (rank == 0) {  
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
      // 
      //  -the minimum distance absmin between 2 particles during the run of the simulation
      //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
      //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
      //
      //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
      //
      printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
      if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
      if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //  
    // Printing summary data
    //  
    if( fsum)
      fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
  }

  //
  //  release resources
  //
  if ( fsum )
    fclose( fsum );
  //free( partition_offsets );
  //free( partition_sizes );
  //free( local );
  free( particles );
  if( fsave )
    fclose( fsave );

  MPI_Finalize( );

  return 0;
}
