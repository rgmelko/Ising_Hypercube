// A program simulate the Ising hamiltonian on a D-dimensional Hypercube
// Roger Melko, June 8, 2013

#include "hypercube.h"
#include "MersenneTwister.h"
#include "simparam.h"
#include <iostream>
using namespace std;

int main(){

    PARAMS param; //read parameter file
    //param.print();

    MTRand mrand(param.SEED_); //random number for metropolis

    HyperCube cube(param.nX_,param.Dim_); //initialize the lattice
    cube.print();

    return 0;

}
