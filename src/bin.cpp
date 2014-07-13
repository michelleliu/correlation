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
#include <assert.h>
#include "common.h"
#include <tuple>

int get_bin_idx(int x_block, int y_block, int z_block, const int mb, const int nb) {
  return (x_block) * (mb) * (nb) + (y_block) * (nb) + (z_block);
}
std::tuple<int, int, int> get_bin_coords(int bin_idx, const int mb, const int nb) {
    int z_block = bin_idx % nb;
    int y_block = ( (bin_idx - z_block) / nb ) % mb;
    int x_block = bin_idx / (mb * nb);
    return std::make_tuple(x_block,y_block,z_block);
}

void get_neighbors(const std::tuple<int, int, int>& block_id, std::vector<std::tuple<int, int, int> >* neighbors) {
    using namespace std;
    printf("For block (%d, %d, %d) nbrs: ", get<0>(block_id), get<1>(block_id), get<2>(block_id));
    for (int i = get<0>(block_id) - 1; i <= get<0>(block_id) + 1; ++i) {
        for (int j = get<1>(block_id) - 1; j <= get<1>(block_id) + 1; ++j) {
            for (int k = get<2>(block_id) - 1; k <= get<2>(block_id) + 1; ++k) {
                if (!(i == get<0>(block_id) && j == get<1>(block_id) && k == get<2>(block_id) )) {
                    if ( i >= 0 && j >= 0 && k >= 0 ) {
                        printf(" (%d, %d, %d), ", i, j, k);
                        neighbors->push_back(make_tuple(i, j, k));
                    }
                }
            }
        }
    }
    printf("\n");
}

int main(int argc, const char * argv[]) {
    using namespace std;
    int mb = 3;
    int nb = 2;
    int bin_idx=0;
    int test_x, test_y, test_z;
    for (int i=0; i < mb; ++i) {
        for (int j=0; j < mb; ++j) {
            for (int k=0; k < nb; ++k) {
                bin_idx=get_bin_idx(i,j,k,mb,nb);
                tuple<int, int, int> bin_coords = get_bin_coords(bin_idx,mb,nb);
                test_x=get<0>(bin_coords);
                test_y=get<1>(bin_coords);
                test_z=get<2>(bin_coords);
                cout << "x, y, z = " << i << ", " << j << ", " << k << endl;
                cout << "bin_idx = " << bin_idx << endl;
                cout << "new x, y, z = " << test_x <<  ", " << test_y << ", " << test_z << endl;
                cout << "--------------------------------" << endl;
            }
        }
    }

    tuple<int, int, int> test_block_id = make_tuple(0, 1, 1);
    vector<tuple<int, int, int> > neighbors;
    get_neighbors(test_block_id,&neighbors);

}
