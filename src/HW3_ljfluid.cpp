//
//  main.cpp
//  HW3_ljfluid
//
//  Created by Michelle Liu on 12/4/13.
//  Copyright (c) 2013 Michelle Liu. All rights reserved.
//

#include <fstream>
#include <iostream>
#include <cmath>
#include <array>
#include <time.h>
#include <stdlib.h>

// Forward declaration
class MDstate;

//  G L O B A L   V A R I A B L E S
//
const int NumAtoms = 108;
std::array < std::array<long double, 3> , NumAtoms> AtomCoords;
double Density = 0.8442;
long double BoxLength = pow((NumAtoms/Density), 1./3.); // Cubic box with the volume sized according to the density
//
//  G L O B A L   V A R I A B L E S

// Reads in coordinates from file
int readCoordinates()
{
    using namespace std;
    ifstream coordinates;
    if (!coordinates.is_open()) {
        coordinates.open("/Users/michelle/Dropbox/BE243/HW3_ljfluid/LJ_108_1.txt");
    }
    if (coordinates.is_open()) {
        for (int iii = 0; iii < NumAtoms; iii++) {
            for (int jjj = 0; jjj < 3; jjj++) {
                coordinates >> AtomCoords[iii][jjj];
            }
        }
    }
    else {
        cout << "Error opening file";
    }
    coordinates.close();
    return 0;
}

// Prints coordinates of all 108 atoms in 2 x 2 array form
int printCoordinates()
{
    using namespace std;
    for (int iii = 0; iii < NumAtoms; iii++) {
        for (int jjj = 0; jjj < 3; jjj++) {
            cout << AtomCoords[iii][jjj] << " ";
        }
        cout << endl;
    }
    return 0;
}

// Distance (scalar) between two particles
// Directly calculated using atom coordinates
long double distance(int first, int second)
{
    long double xDist = AtomCoords[second][0] - AtomCoords[first][0];
    long double yDist = AtomCoords[second][1] - AtomCoords[first][1];
    long double zDist = AtomCoords[second][2] - AtomCoords[first][2];
    if (xDist > BoxLength/2) {
        xDist = BoxLength - xDist;
    }
    if (yDist > BoxLength/2) {
        yDist = BoxLength - yDist;
    }
    if (zDist > BoxLength/2) {
        zDist = BoxLength - zDist;
    }
    long double Dist = sqrtl(xDist*xDist + yDist*yDist + zDist*zDist);
    return Dist;
}

// Potential between two particles
long double pairPotential(int first, int second)
{
    long double Dist = distance(first, second);
    long double Vr = 4*(pow((1/Dist),12) - pow((1/Dist),6));
    return Vr;
}

// Total potential on a single particle from all other particles
long double particlePotential(int particleIndex)
{
    long double Potential = 0;
    for (int iii = 0; iii < NumAtoms; iii++) {
        if (iii != particleIndex) {
            Potential = Potential + pairPotential(particleIndex, iii);
        }
    }
    return Potential;
}

// Total potential energy of atom cluster
long double totalPotential()
{
    long double TotPotential = 0;
    for (int iii = 0; iii < NumAtoms; iii++) {
        for (int jjj = iii+1; jjj < NumAtoms; jjj++) {
            TotPotential = TotPotential + pairPotential(iii, jjj);
        }
    }
    return TotPotential;
}

// Vector pointing to first atom from the second atom
std::array<long double, 3> pairVector(int first, int second)
{
    long double xDist = AtomCoords[second][0] - AtomCoords[first][0];
    long double yDist = AtomCoords[second][1] - AtomCoords[first][1];
    long double zDist = AtomCoords[second][2] - AtomCoords[first][2];
    if (xDist > BoxLength/2) {
        xDist = BoxLength - xDist;
    }
    if (yDist > BoxLength/2) {
        yDist = BoxLength - yDist;
    }
    if (zDist > BoxLength/2) {
        zDist = BoxLength - zDist;
    }
    std::array<long double, 3> Vector = {xDist, yDist, zDist};
    return Vector;
}

// Length of vector pointing to first atom from second atom (same thing as distance function)
long double vectorLength(std::array<long double, 3> Vector)
{
    long double Dist = sqrt(Vector[0]*Vector[0] + Vector[1]*Vector[1] + Vector[2]*Vector[2]);
    return Dist;
}

// Vector force on first atom from second atom
std::array<long double, 3> pairForce(int first, int second)
{
    std::array<long double, 3> PairForce = {0, 0, 0};
    std::array<long double, 3> Vector = pairVector(first, second);
    long double Distance = vectorLength(Vector);
    long double ForceMagnitude;
    ForceMagnitude = 4*(-12*pow(1/Distance, 13) + 6*pow(1/Distance, 7));
    for (int iii = 0; iii < 3; iii++) {
        PairForce[iii] = (ForceMagnitude*Vector[iii])/Distance;
    }
    return PairForce;
}

// Total vector force on one atom from all other atoms
std::array<long double, 3> particleForce(int particleIndex)
{
    std::array<long double, 3> TotalForce = {0,0,0};
    std::array<long double, 3> PairForce = {0,0,0};
    for (int iii = 0; iii < NumAtoms; iii++) {
        if (iii != particleIndex) {
            PairForce = pairForce(particleIndex, iii);
            for (int jjj = 0; jjj < 3; jjj++) {
                TotalForce[jjj] = TotalForce[jjj] + PairForce[jjj];
            }
        }
    }
    return TotalForce;
}

// Generate three uniform random numbers
// Generate three M-B distro random numbers
// Use to initialize velocities (only used at the beginning)
std::array<long double,3> randVelocities()
{
    std::array<long double, 3> Velocities = {0,0,0};
    std::array<long double, 4> Uniform = {0,0,0,0};
    long double x1, x2, r, x3, x4, r2;
    // initialize random seed
    srand(time(NULL));
    
    do {
        for (int iii = 0; iii < 2; iii++) {
            Uniform[iii] = 1.0*rand()/RAND_MAX;
            //std::cout << Uniform[iii] << std::endl;
        }
        r = Uniform[0]*Uniform[0] + Uniform[1]*Uniform[1];
    } while (r > 1.0);
    
    do {
        for (int iii = 2; iii < 4; iii++) {
            Uniform[iii] = 1.0*rand()/RAND_MAX;
            //std::cout << Uniform[iii] << std::endl;
        }
        r2 = Uniform[2]*Uniform[2] + Uniform[3]*Uniform[3];
    } while (r2 > 1.0);
    
    x1 = 2.0*Uniform[0] - 1.0;
    x2 = 2.0*Uniform[1] - 1.0;
    x3 = 2.0*Uniform[2] - 1.0;
    x4 = 2.0*Uniform[3] - 1.0;
    
    /*
    std::cout << "r = " << r << std::endl;
    std::cout << "r2 = " << r2 << std::endl;
    */
    
    r = sqrtl((-2.0*logl(r))/r);
    r2 = sqrtl((-2.0*logl(r2))/r2);
    Velocities[0] = x1*r;
    Velocities[1] = x2*r;
    Velocities[2] = x3*r2;
    
    
    return Velocities;
}

// Create a class and name it "State"
class MDstate
{
public:
    double tstep;
    std::array < std::array<long double, 3> , NumAtoms> Position;
    std::array < std::array<long double, 3> , NumAtoms> Velocity;
    std::array < std::array<long double, 3> , NumAtoms> Force;
    int Ntimesteps = 0;
    
    void MDinitialize(double delt, std::array < std::array<long double, 3> , NumAtoms> pos, std::array < std::array<long double, 3> , NumAtoms> vel, std::array < std::array<long double, 3> , NumAtoms> forces, int Nsteps)
    {
        tstep = delt;
        Position = pos;
        Velocity = vel;
        Force = forces;
        Ntimesteps = Nsteps;
    }
    
    void MDmove()
    {
        using namespace std;
        if (Ntimesteps == 0) {
            cout << "MD state not initialized" << endl;
            exit(0);
        }
        else
        {
            // Velocity Verlet Algorithm
            
        }
    }
};


/* I don't think I need to do this
// Potential correction
// What is this even for?
long double potentialCorrection(long double Potential)
{
    long double CorrectedPotential;
    long double Rc;
    return 0;
}
*/

MDstate SystemState;

int main(int argc, const char * argv[])
{
    using namespace std;
    readCoordinates();
    cout << "Box length = " << BoxLength << endl;
    cout << "Potential between particles 29 and 30 = " << pairPotential(29, 30) << endl;
    cout << "Total potential on particle 30 = " << particlePotential(30) << endl;
    
    // Calculate force on particle 30
    std::array<long double, 3> TotalForce = particleForce(30);
    cout << "Total (vector) force on particle 30 = " << endl;
    cout << "     ";
    for (int iii = 0; iii < 3; iii ++) {
        cout << TotalForce[iii] << " " ;
    }
    cout << endl;

    randVelocities();
    
    SystemState.MDmove();
    // printCoordinates();
    return 0;
}