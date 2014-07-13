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

int get_bin_idx(int x_block, int y_block, int z_block, const int mb, const int nb) {
  return (x_block) * (mb) * (nb) + (y_block) * (nb) + (z_block);
}
std::tuple<int, int, int> get_row_col(int bin_idx, const int mb, const int nb) {
    int x_block = bin_idx % (mb * nb);
    int y_block = bin_idx;
    int z_block = bin_idx % nb;
    return std::make_tuple(x_block,y_block,z_block);
}

int main(int argc, const char * argv[]) {
    int mb = 3;
    int nb = 2;
    int bin_idx=0;
    int test_x, test_y, test_z;
    for (int i=0; i < mb; ++i) {
        for (int j=0; j < mb; ++j) {
            for (int k=0; k < nb; ++k) {
                bin_idx=get_bin_idx(i,j,k,mb,nb);
                std::tuple<int, int, int> bin_coords = get_bin_coords(bin_idx,mb,nb);
                test_x=std::get<0>(bin_coords);
                test_y=std::get<1>(bin_coords);
                test_z=std::get<2>(bin_coords);
                std::cout << "x, y, z = " << i << ", " << j << ", " << k << endl;
                std::cout << "bin_idx = " << bin_idx << endl;
                std::cout << "new x, y, z = " << test_x <<  ", " << test_y << ", " << test_z << endl;
            }
        }
    }

}
