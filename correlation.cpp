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

int main(int argc, const char * argv[]) {
    using namespace std;
    filebuf fb;
    if (fb.open("example22.traj",std::ios::in))
    {
        istream is (&fb);
        while (is)
            cout << char(is.get());
        fb.close();
    }
    return 0;
}
