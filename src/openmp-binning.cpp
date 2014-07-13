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

// typedef unordered_map<pair<int, int>, vector<int>*, pairhash> particle_map;
typedef unordered_map<pair<int, int>, vector<int>*, pairhash>::iterator particle_map_itr;

bool is_padding_bin(pair<int, int> block_id, int num_bins) {
  return !(block_id.first >= 0 && block_id.first < num_bins &&
          block_id.second >= 0 && block_id.second < num_bins);
}

// void print_bins(particle_map* bins, particle_t* particles, int num_bins) {
//   for (particle_map_itr itr = bins->begin(); itr != bins->end(); ++itr) {
//     if (!is_padding_bin(itr->first, num_bins)) {
//       printf("Bin (%d, %d): ", itr->first.first, itr->first.second);
//       for (int i = 0; i < itr->second->size(); ++i) {
//         printf("%d, ", itr->second->at(i));
//       }
//       printf("\n");
//     }
//   }
// }

// void clear_bins(particle_map* bins) {
//   for (particle_map_itr itr = bins->begin(); itr != bins->end(); ++itr) {
//     itr->second->clear();
//   }
// }

// Populate bins from list of particles.
// Puts each particle into <block_row, block_col> bin
// void populate_bins(particle_map* bins, particle_t* particles, int num_particles) {
//   for (int i = 0; i < num_particles; ++i) {
//     // Get x, y co-ordinates of this particle
//     int row_block = (int) floor(particles[i].x / (BIN_SIZE));
//     int col_block = (int) floor(particles[i].y / (BIN_SIZE));
//     pair<int, int> block_id(row_block, col_block);
//     particle_map_itr itr = bins->find(block_id);
//     if (itr != bins->end()) {
//       itr->second->push_back(i);
//     }
//   }
// }

int get_bin_idx(int bin_row, int bin_col, const int num_bins) {
  return (bin_row + 1) * (num_bins + 2) + (bin_col + 1);
}

pair<int, int> get_row_col(int bin_idx, const int num_bins) {
  int bin_col = bin_idx % (num_bins + 2) - 1;
  int bin_row = bin_idx / (num_bins + 2) - 1;
  return make_pair(bin_row, bin_col);
}

void populate_bins(vector<int>** bins, particle_t* particles, int num_particles, omp_lock_t* locks, int num_bins) {
  #pragma omp for
  for (int i = 0; i < num_particles; ++i) {
    // Get x, y co-ordinates of this particle
    int row_block = (int) floor(particles[i].x / (BIN_SIZE));
    int col_block = (int) floor(particles[i].y / (BIN_SIZE));

    int id = get_bin_idx(row_block, col_block, num_bins);
    omp_set_lock(&locks[id]);
    bins[id]->push_back(i);
    omp_unset_lock(&locks[id]);

    // pair<int, int> block_id(row_block, col_block);
    // particle_map_itr itr = bins->find(block_id);
    // if (itr != bins->end()) {
    //   itr->second->push_back(i);
    // }
  }
}

void clear_bins(vector<int>** bins, omp_lock_t* locks, const int num_bins) {
  for (int i = 0; i < (num_bins + 2)*(num_bins+2); ++i) {
    //omp_set_lock(&locks[i]);
    bins[i]->clear();
    //omp_unset_lock(&locks[i]);
  }
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
  // printf("\n");
}

void apply_force_bins(vector<int>* source, vector<int>* other, particle_t* particles,
    double *dmin, double *davg, int *navg) {
  for (int iaf = 0; iaf < source->size(); ++iaf) {
    for (int jaf = 0; jaf < other->size(); ++jaf) {
      // printf("Applying force to %d from %d\n", source->at(i), other->at(j));
      apply_force(particles[source->at(iaf)],
                  particles[other->at(jaf)],
                  dmin, davg, navg);
    }
  }
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{
  int navg,nabsavg=0, numthreads;
  double davg,dmin, absmin=1.0, absavg=0.0;

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

  FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
  FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

  particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
  set_size( n );
  init_particles( n, particles );

  // Create a hash_map from bin_location to set of particles in a bin
  // particle_map bins;
  const int num_bins = ceil(get_size() / (BIN_SIZE));

  vector<int>** bins;
  bins = new vector<int>*[(num_bins+2)*(num_bins+2)];
  omp_lock_t* locks;
  locks = new omp_lock_t[(num_bins+2)*(num_bins+2)];

  // vector<pair<int, int> > bin_idxs;
  printf("Size is %lf, bin size is %lf, Number of bins %d\n", get_size(), BIN_SIZE,
      num_bins);

  for (int i = -1; i < num_bins + 1; ++i) {
    for (int j = -1; j < num_bins + 1; ++j) {
      // pair<int, int> mypair = make_pair(i,j);
      // bin_idxs.push_back(mypair);
      int id = get_bin_idx(i, j, num_bins);
      bins[id] = new vector<int>();
      omp_init_lock(&locks[id]);
    }
  }

  //
  //  simulate a number of time steps
  //
  double simulation_time = read_timer( );

  // Count num particles in each bin
  // print_bins(&bins, particles, num_bins);

  #pragma omp parallel shared(bins,locks) private(dmin)
  {
  numthreads = omp_get_num_threads();

  populate_bins(bins, particles, n, locks, num_bins);
  for( int step = 0; step < NSTEPS; step++ )
  {
    navg = 0;
    davg = 0.0;
    dmin = 1.0;

    //
    //  compute forces
    //

    // for bin in bins
    //    get neigboring bins
    //    for each neigboring bin
    //      do apply_force
    #pragma omp for reduction (+:navg) reduction(+:davg)
    for (int j = 0; j< (num_bins+2)*(num_bins+2); ++j) {
      // particle_map_itr itr = bins.find(bin_idxs[j]);
      // vector<int>* this_bin = itr->second;
      vector<int>* this_bin = bins[j];
      // List surrounding bins using key
      pair<int, int> block_id = get_row_col(j, num_bins);

      if (!is_padding_bin(block_id, num_bins)) {
        for (int i = 0; i < this_bin->size(); ++i) {
          particles[this_bin->at(i)].ax = particles[this_bin->at(i)].ay = 0;
        }

        apply_force_bins(this_bin, this_bin, particles, &dmin, &davg, &navg);
        vector<pair<int, int> > neighbors;
        get_neighbors(block_id, &neighbors);

        for (vector<pair<int, int> >::iterator nbr_itr = neighbors.begin();
             nbr_itr != neighbors.end(); ++nbr_itr) {
          vector<int>* nbr_bin = bins[get_bin_idx(nbr_itr->first, nbr_itr->second, num_bins)];
          apply_force_bins(this_bin, nbr_bin, particles, &dmin, &davg, &navg);
        }
      }
    }

    //
    //  move particles
    //
    #pragma omp for
    for( int i = 0; i < n; i++ ){
      move( particles[i] );
    }

    // Rebin particles here
    #pragma omp master
    {
    clear_bins(bins, locks, num_bins);
    }
    #pragma omp barrier
    populate_bins(bins, particles, n, locks, num_bins);
    // printf("After step %d\n", step);
    // Count num particles in each bin
    // print_bins(&bins, particles, num_bins);


    if( find_option( argc, argv, "-no" ) == -1 )
    {
      //
      // Computing statistical data
      //
      #pragma omp master
      if (navg) {
        absavg +=  davg/navg;
        nabsavg++;
      }
      #pragma omp critical
      if (dmin < absmin) absmin = dmin;

      //
      //  save if necessary
      //
      #pragma omp master
      if( fsave && (step%SAVEFREQ) == 0 )
        save( fsave, n, particles );
    }

    // NOTE: Barrier here to make sure all threads are in the same step
    #pragma omp barrier
  }
  }
  simulation_time = read_timer( ) - simulation_time;

  printf( "threads %d, n = %d, simulation time = %g seconds", numthreads, n, simulation_time);

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
    fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

  //
  // Clearing space
  //
  if( fsum )
    fclose( fsum );
  free( particles );
  if( fsave )
    fclose( fsave );

  return 0;
}
